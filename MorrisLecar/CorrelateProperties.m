function resultS = CorrelateProperties(xData, yData, options, varargin)

defaultOptions = {'fitPooled', true, 'fitSubTypes', false, ...
		  'compareDist', false};
%Note:  I don't remember what compareDist does!!!
options = GetOptions(defaultOptions, options, true);

%Setup resultS
if length(varargin) >= 1
  shuffleTrial = true;
  
  resultS = varargin{1};
  xData = resultS.xData;
  yData = resultS.yData;
else
  shuffleTrial = false;
  
  resultS = constructResultS(xData, yData, options);
end
resultS.current = 1;

%First do the correlation with all cells grouped together
fitNum = 1;
if options.fitPooled
  label = 'Pooled';
  useInd = 1:size(xData.mat, 1);  % use all cells
  resultS = mainCorrIPRoutine(options, resultS, shuffleTrial, ...
			      label, fitNum);
  fitNum = fitNum + 1;
end

%Now separate the cells out by type:
if options.fitSubTypes
  for n = 1:length(resultS.cellTypes)
    cellType = resultS.cellTypes{n};
    label = cellType;
    useInd = find(strcmp(xData.cellType, cellType));
    resultS = mainCorrIPRoutine(options, resultS, shuffleTrial, ...
				label, fitNum);
    fitNum = fitNum + 1;
  end
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xData, yData] = removeCrap(xData, yData)
%Remove all the rows where XData is all NaN or yData is all NaN
%  (actually, keep the ones where that's not true)
keepRows = find( ...
    sum(isfinite(xData.mat), 2) > 0 & ...
    sum(isfinite(yData.mat), 2) > 0 ...
    );

xData.mat = xData.mat(keepRows,:);
xData.cellType = xData.cellType(keepRows);
xData.ID = xData.ID(keepRows);
yData.mat = yData.mat(keepRows,:);
yData.cellType = yData.cellType(keepRows);
yData.ID = yData.ID(keepRows);
%fprintf('Size of %s: %gx%g\n', xData.label, size(xData.mat, 1), ...
%	size(xData.mat, 2))
%fprintf('Size of %s: %gx%g\n', yData.label, size(yData.mat, 1), ...
%	size(yData.mat, 2))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resultS = constructResultS(xData, yData, options)
%Remove all the rows where XData is all NaN or yData is all NaN
[xData, yData] = removeCrap(xData, yData);
resultS.label = sprintf('%s vs %s', xData.label, yData.label);
resultS.xData = xData;
resultS.yData = yData;

%Make list all the cell types:
resultS.cellTypes = unique({resultS.xData.cellType{:}, ...
		    resultS.yData.cellType{:}});
numCellTypes = length(resultS.cellTypes);

numTypeFits = options.fitPooled + ...
    numCellTypes * options.fitSubTypes;
numX = size(xData.mat, 2);  %this is the number of x properties,
			    %not cells
numY = size(yData.mat, 2);  %this is the number of y properties

%preliminary estimate of numFits.
numFits = numTypeFits * (numX + numY);

useRows = cell(numTypeFits, 1);
useXCols = cell(numTypeFits, 1);
useYCols = cell(numTypeFits, 1);
numFits = 0;
fitNum = 1;
if options.fitPooled
  useRows{fitNum} = 1:size(xData.mat);
  useXCols{fitNum} = find(sum(isfinite(xData.mat), 1) > 0);
  useYCols{fitNum} = find(sum(isfinite(yData.mat), 1) > 0);
  numFits = numFits + length(useXCols{fitNum}) + ...
	    length(useYCols{fitNum});
  fitNum = fitNum + 1;
end
if options.fitSubTypes
  for n = 1:numCellTypes
    cellType = resultS.cellTypes{n};
    useRows{fitNum} = find(strcmp(xData.cellType, cellType));
    useXCols{fitNum} = ...
	find(sum(isfinite(xData.mat(useRows{fitNum},:)), 1) > 0);
    useYCols{fitNum} = ...
	find(sum(isfinite(yData.mat(useRows{fitNum},:)), 1) > 0);
    numFits = numFits + length(useXCols{fitNum}) + ...
	      length(useYCols{fitNum});
    fitNum = fitNum + 1;
  end
end
resultS.useRows = useRows;
resultS.useXCols = useXCols;
resultS.useYCols = useYCols;

resultS.numFits = numFits;
resultS.fitErr = zeros(numFits, 1);
resultS.RSquared = zeros(numFits, 1);
resultS.descriptions = cell(numFits, 1);
resultS.indices = cell(numFits, 1);
resultS.y = cell(numFits, 1);
resultS.yPredict = cell(numFits, 1);
resultS.coefs = cell(numFits, 1);
resultS.pVals = zeros(numFits, 1);
resultS.sections = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resultS = mainCorrIPRoutine(options, resultS, shuffleTrial, ...
				     label, fitNum)

useRows = resultS.useRows{fitNum};
xCols = resultS.useXCols{fitNum};
yCols = resultS.useYCols{fitNum};
numCells = length(useRows);
numXTypes = length(xCols);
numYTypes = length(yCols);

if(shuffleTrial)
  xRows = useRows(randperm(numCells));
  x = resultS.xData.mat(xRows, xCols);
  yRows = useRows(randperm(numCells));
  y = resultS.yData.mat(yRows, yCols);
  
  shuffleNum = resultS.shuffleNum;
else
  x = resultS.xData.mat(useRows, xCols);
  y = resultS.yData.mat(useRows, yCols);
  
  xDataLabel = resultS.xData.label;
  yDataLabel = resultS.yData.label;
  xLabels = cell(numXTypes, 1);
  yLabels = cell(numYTypes, 1);
  for n = 1:numXTypes
    xLabels{n} = [label, ': ', resultS.xData.labels{xCols(n)}];
  end
  for n = 1:numYTypes
    yLabels{n} = [label, ': ', resultS.yData.labels{yCols(n)}];
  end

  %What to do with this?
  %if(options.IPvsIP)
  %  IPvsIP(X, XLabels, Label, scoreHandle, options)
  %end
  if(options.compareDist)
    CompareXDistYDist(x, y, Label);
  end
end

%Fit yData with linear combinations of xData
n1 = resultS.current;
n2 = n1 + numYTypes - 1;
if(shuffleTrial)
  resultS.shuffleFitErr(n1:n2,shuffleNum) = FindCorrelations(x, y);
else
  [resultS.fitErr(n1:n2), resultS.RSquared(n1:n2), ...
   resultS.indices(n1:n2), resultS.y(n1:n2), ...
   resultS.yPredict(n1:n2), resultS.coefs(n1:n2), ...
   resultS.pVals(n1:n2)] = FindCorrelations(x, y);
  resultS.descriptions(n1:n2) = getDescriptions(xLabels, yLabels, ...
						resultS.indices(n1:n2));
  section.label = sprintf('%s: Fitting %s with %s', ...
			  label, yDataLabel, xDataLabel);
  section.inds = n1:n2;
  section.fitLabels = yLabels;
  section.paramLabels = stripCellID(xLabels);
  resultS.sections = [resultS.sections, section];
end

%Fit xData with linear combinations of xData
n1 = n2 + 1;
n2 = n2 + numXTypes;
if(shuffleTrial)
  resultS.shuffleFitErr(n1:n2,shuffleNum) = FindCorrelations(y, x);
else
  [resultS.fitErr(n1:n2), resultS.RSquared(n1:n2), ...
   resultS.indices(n1:n2), resultS.y(n1:n2), ...
   resultS.yPredict(n1:n2), resultS.coefs(n1:n2), ...
   resultS.pVals(n1:n2)] = FindCorrelations(y, x);
  resultS.descriptions(n1:n2) = getDescriptions(yLabels, xLabels, ...
						resultS.indices(n1:n2));
  section.label = sprintf('%s: Fitting %s with %s', ...
			  label, xDataLabel, yDataLabel);
  section.inds = n1:n2;
  section.fitLabels = xLabels;
  section.paramLabels = stripCellID(yLabels);
  resultS.sections = [resultS.sections, section];
end

resultS.current = n2 + 1;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function descriptions = getDescriptions(paramLabels, fitLabels, indices)
paramLabels = stripCellID(paramLabels);  %Don't want to keep
                                         %repeating cell ID

numFit = length(fitLabels);
descriptions = cell(numFit, 1);
for n = 1:numFit
  descriptions{n} = sprintf('%s fit with', fitLabels{n});
  for ind = indices{n}
    descriptions{n} = sprintf('%s %s', descriptions{n}, ...
			      paramLabels{ind});
  end
end
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Labels = stripCellID(Labels)
if(strcmp(class(Labels), 'cell'))
  for n = 1:length(Labels)
    Temp = Labels{n};
    Ind = strfind(Temp, ':');
    if(length(Ind) > 0)
      Temp = Temp((Ind+1):end);
      Labels{n} = Temp;
    end
  end
elseif(strcmp(class(Labels), 'char'))
  Ind = strfind(Labels, ':');
  if(length(Ind) > 0)
    Labels = Labels((Ind+2):end);
  end
else
  error('Labels must be string or cell array of strings')
end
return

%%%%%%%%% The following should be moved elsewhere %%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IPvsIP(X, XLabels, Label, scoreHandle, options)
TitleBase = [Label, ' IP'];
Title = [Label, ' IP Correlations'];

fprintf('%s vs IP p-values:', TitleBase)
NumX = size(X, 2);
RSqrd = zeros(NumX, 1);
numFits = NumX * NumX * (options.fitPooled + options.fitSubTypes);

for n = 1:NumX
  X_n = X(:,n);
  XOther = X(:,[1:(n-1), (n+1):end]);
  [XFit, XPVal] = JackknifeFit(XOther, X_n);
  
  Err = XFit - X_n;
  Var = X_n - mean(X_n);
  RSqrd(n) = 1 - Err' * Err / (Var' * Var);
  fprintf(' %s=%.2g', stripCellID(XLabels{n}), XPVal * numFits)
end
fprintf('\n')
fprintf('%s vs IP R^2:', TitleBase)
for n = 1:NumX
  fprintf(' %s=%.2g', stripCellID(XLabels{n}), RSqrd(n))
end
fprintf('\n')

RSqrdMat = zeros(NumX, NumX);
PMat = zeros(NumX, NumX);
fprintf('%s Pairwise R^2:\n', TitleBase)
for n = 1:NumX
  for m = 1:NumX
    X_n = X(:,n) - mean(X(:,n));
    X_m = X(:,m) - mean(X(:,m));

    [RSqrdMat(n,m), PMat(n,m)] = GetLineFit(X_n, X_m);
    fprintf('\t%.2g', RSqrdMat(n,m))
  end
  fprintf('\n')
end

fprintf('%s Pairwise p-values:\n', TitleBase)
for n = 1:NumX
  for m = 1:NumX
    fprintf('\t%.2g', PMat(n,m) * numFits)
  end
  fprintf('\n')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RSqrd, PVal] = GetLineFit(X_n, X_m)
Ind = find(isfinite(X_n) & isfinite(X_m));
X_n = X_n(Ind);
X_m = X_m(Ind);
Coefs = polyfit(X_n, X_m, 1);
Fit_m = Coefs(1) * X_n + Coefs(2);

Err = Fit_m - X_m;
Var = X_m - mean(X_m);
RSqrd = 1 - Err' * Err / (Var' * Var);
Fit_m = JackknifeFit(X_n, X_m);
[PearsonR, PVal] = Pearson(X_m, Fit_m);
return
