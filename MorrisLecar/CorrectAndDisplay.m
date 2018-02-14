function CorrectAndDisplay(resultS, options)

if canShuffleCorrect(resultS)
  resultS = shuffleCorrect(resultS);
elseif isfield(resultS, 'pVals')
  resultS = holmBonferroniCorrect(resultS);
else
  error(['Can''t correct because there are no p-values or', ...
	 ' shuffled data.\n']);
end

displayCorrectedResults(resultS, options);

displaySummary(resultS);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boolVal = canShuffleCorrect(resultS)
boolVal = isfield(resultS, 'shuffleFitErr') && ...
	  isfield(resultS, 'fitErr') && ...
	  size(resultS.shuffleFitErr, 2) > 0;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resultS = shuffleCorrect(resultS, sigLevel)
if(nargin < 2)
  sigLevel = 0.05;
end

numShuffle = size(resultS.shuffleFitErr, 2);
numFits = resultS.numFits;
resultS.sigLevels = zeros(numFits, 1); %stores max values of fitErr that
				       %  are significant
resultS.sig = repmat(false, numFits, 1); %is significan true/false for each
					 %  fit
resultS.pAdjusted = ones(numFits, 1);  %stores p-values adjusted for
				       %multiple comparisons

%Get dataMat, a numFits x numShuffle matrix
dataMat = resultS.shuffleFitErr;
%Sort the fit error from smallest to largest
[sortData, sortInds] = sort(resultS.fitErr);

includeInds = 1:numFits;
for n = 1:numFits
  %Loop through data, starting with the best fit.  Get the best
  %  fits to shuffled data, and see how many are better
  bestFits = min(dataMat(includeInds,:), [], 1);
  ind = sortInds(n);
  numBetter = sum(bestFits < sortData(n));
  %Get a p-value on this fit being better than fits to shuffled data
  resultS.pAdjusted(ind) = (numBetter + 1) / (numShuffle + 1);
  %Get a fit error that would meet the significance criterion
  bestFits = sort(bestFits);
  resultS.sigLevels(ind) = bestFits(round(sigLevel * numShuffle));
  
  if(resultS.pAdjusted(ind) < sigLevel)
    %if the adjusted p-value is small enough, the fit is
    % significant.
    resultS.sig(ind) = true;
    %Remove the shuffled trials of *this* fit from consideration,
    %  and go on to the next-best fit.
    includeInds = setdiff(includeInds, ind);
  else
    %The adjusted p-value is NOT small enough so this fit and all
    %  subsequent fits are NOT significant.
    str = sprintf('Fit %g of %g not significant (p=%g)', ...
		  n, numFits, resultS.pAdjusted(ind));
    fprintf('%s\n', str)
    h = NamedFigure(str);
    set(h, 'WindowStyle', 'docked')
    
    [histNums, histX] = hist(log(bestFits), 100);
    val = log(sortData(n));
    bar(histX, histNums)
    hold on
    plot([val, val], [0, max(histNums)], 'r-')
    hold off
    
    resultS.sigLevels(sortInds((n+1):end)) = resultS.sigLevels(ind);
    resultS.pAdjusted(sortInds((n+1):end)) = 1.0;
    break
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resultS = holmBonferroniCorrect(resultS, sigLevel)
if(nargin < 2)
  sigLevel = 0.05;
end

[sortP, sortInd] = sort(resultS.pVals, 'descend');

numFits = resultS.numFits;
resultS.sigLevels = zeros(numFits, 1);
resultS.sig = repmat(false, numFits, 1);
resultS.pAdjusted = ones(numFits, 1);

for n = numFits:-1:1
  ind_n = sortInd(n);
  resultS.sigLevels(ind_n) = sigLevel / n;
  if(sortP(n) < sigLevels(ind_n))
    resultS.sig(ind_n) = true;
    resultS.pAdjusted(ind_n) = sortP(n) * n;
  else
    resultS.sigLevels(sortInd((n-1):-1:1)) = sigLevels(ind_n);
    break
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayCorrectedResults(resultS, options, sigLevel)
if(nargin < 3)
  sigLevel = 0.05;
end
labelSize = 18;
titleSize = 18;

numSections = length(resultS.sections);
for n = 1:numSections
  section = resultS.sections(n);
  fprintf('%s %s %s\n', '########', section.label, '########')
  for ind = section.inds
    %Report significance
    if(resultS.sig(ind))
      fprintf('%s ', '[SIGNIFICANT]')
    else
      fprintf('%s ', '[fudge]      ')
    end
    
    %Report fit names and R^2/p-values
    fprintf('%s\n', resultS.descriptions{ind})
    fprintf('               p=%g', resultS.pAdjusted(ind))
    if isfield(resultS, 'RSquared')
      fprintf(' R^2=%g', resultS.RSquared(ind))
    end
    if isfield(resultS, 'pVals')
      fprintf(' pRaw=%g', resultS.pVals(ind))
    end
    if isfield(resultS, 'fitErr')
      fprintf(' rChiSquared=%g', resultS.fitErr(ind))
    end
    %resultS.sigLevels?
    fprintf('\n')

    %Make individual plots
    if(options.plotIndividual)
      y = resultS.y{ind};
      yPredict = resultS.yPredict{ind};
      yRange = [min(min(y), min(yPredict)), ...
		max(max(y), max(yPredict))];
      titleStr = resultS.descriptions{ind};
      yLabelInd = strfind(titleStr, ' fit');
      yLabel = titleStr(1:(yLabelInd-1));
      h = NamedFigure(titleStr);
      set(h, 'WindowStyle', 'docked')
      hold off
      plot(yRange, yRange, 'r-')
      hold on
      plot(y, yPredict, 'b.', 'MarkerSize', 20)
      hold off
      xlabel(['Scored ', yLabel], 'FontName', 'Arial', ...
	     'FontSize', labelSize);
      ylabel(['Fit to Scored ', yLabel], 'FontName', 'Arial', ...
	     'FontSize', labelSize);
      title(titleStr, 'FontName', 'Arial', 'FontSize', titleSize);
    end
  end
  %Make overall section plots
  if(options.plotImportance)
    numY = length(section.fitLabels);
    numX = length(section.paramLabels);
    displayMat = zeros(numY, numX);
    for m = 1:numY
      ind = section.inds(m);
      displayMat(m, resultS.indices{ind}) = -log10(resultS.pAdjusted(ind));
    end
    titleStr = section.label;
    h = NamedFigure(titleStr);
    set(h, 'WindowStyle', 'docked')
    FigureImportance(displayMat, section.fitLabels, ...
		     section.paramLabels, 0.10);
    title(titleStr, 'FontName', 'Arial', 'FontSize', titleSize);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displaySummary(resultS, sigLevel)
if(nargin < 2)
  sigLevel = 0.05;
end

labels = getAllLabels(resultS);
numLabels = length(labels);
numSig = zeros(numLabels, 1);

for n = 1:length(resultS.descriptions)
  if(resultS.sig(n))
    desc = resultS.descriptions{n};
    for m=1:numLabels
      if(StringCheck(desc, labels{m}))
	numSig(m) = numSig(m) + 1;
      end
    end
  end
end

[numSig, sortInd] = sort(numSig, 'descend');
labels = labels(sortInd);

disp('Number of significant fits involving each property:')
for n = 1:length(labels)
  fprintf('  %s:  %g\n', labels{n}, numSig(n))
end

pVal = resultS.pAdjusted;
n = histc(pVal, [0, 0.001, 0.01, 0.05]);
n(4) = length(pVal) - sum(n(1:3));
n = n / length(pVal);

titleStr = ['Distribution of adjusted p-values for ', resultS.label];
h = NamedFigure(titleStr);
set(h, 'WindowStyle', 'docked')
bar(n)
set(gca, 'XTickLabel', {'< .001', ' < .01', ' < .05', '>= .05'})
ylabel('Fraction of fits', 'FontName', 'Arial', 'FontSize', 18)
xlabel('p-values', 'FontName', 'Arial', 'FontSize', 18)
title(titleStr, 'FontName', 'Arial', 'FontSize', 18)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels = getAllLabels(resultS)
labels = {};
for n = 1:length(resultS.descriptions)
  desc = resultS.descriptions{n};
  n1 = strfind(desc, ': ');
  if(length(n1) == 0)
    n1 = 1;
  else
    n1 = n1(1) + 2;
  end
  n2 = strfind(desc, 'fit with');
  n2 = n2(1) - 1;
  label = desc(n1:n2);
  
  labels = {labels{:}, label};
end

labels = unique(labels);
return
