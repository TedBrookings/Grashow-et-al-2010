function [rChiSquared, RSquared, varargout] = ...
    FindCorrelations(x, y, varargin)
% [rChiSquared, RSquared, fitIndices, yCells, yPredict, ...
%  coefs, pVals] = FindCorrelations(x, y, maxLookAhead)
%
%  Try to find correlations of the form y = A * x + B
%  Search through relationships with varying numbers of xs,
%  using f-statistic to justify adding more.
%  The first form will print out values or make plots, while the
%  second form will silently return the results.
%
%  INPUTS:
%    x: Dependent variables, numData x numX matrix
%    y: Independent variables, numData x numY matrix
%   (OPTIONAL)
%     maxLookAhead: number of parameters to consider adding at once
%                   to the current best fit.  Defaults to 3.
%     maxCombos: maximum number of parameters to fit with.
%                Defaults to Inf
%  OUTPUTS:
%    rChiSquared: numY x 1 array of reduced sum of chi^2 per degree
%                 of freedom for each y fit
%    RSquared:  numY x 1 array of R^2 for each y fit
%   (OPTIONAL)
%    fitIndices: numY x 1 cell array of arrays of x-indices.  The
%                nth cell array gives the x-indices used to fit
%                the nth Y variable.
%    yCells: numY x 1 cell array of individial y traces
%    YPredict: numY x 1 cell array of predicted y traces
%    coefs:   numY x 1 cell array of coefficients from fit
%    pVals:  numY x 1 cell array of p-values obtained by using
%            JackknifeFit and a Pearson correlation test.

defaultOptions = {'maxLookAhead', 3, 'maxCombos', Inf};
options = GetOptions(defaultOptions, varargin, true);
if nargout >= 7
  options.getPVals = true;
else
  options.getPVals = false;
end

numData = size(x, 1);
if(size(y,1) ~= numData)
  error('x and y must have the same number of rows')
end
numX = size(x,2);
numY = size(y,2);

rChiSquared = zeros(numY, 1);
RSquared = zeros(numY, 1);
fitIndices = cell(numY, 1);
yCell = cell(numY, 1);
yPredictCell = cell(numY, 1);
coefs = cell(numY, 1);
if options.getPVals
  pVals = zeros(numY, 1);
end

for m = 1:numY
  y_m = y(:,m);
  yCell{m} = y_m;
  AllInd = 1:numX;

  bestInd = [];
  rChiSquared(m) = Inf;
  Done = false;
  numAdd = 1;
  while(~Done)
    %Put the unused indices in Missing
    Missing = setdiff(AllInd, bestInd);
    %Figure out all the combinations of missing indices to try adding:
    if(numAdd == 1)
      Combos = Missing';
    else
      Combos = nchoosek(Missing, numAdd);
    end
    numCombos = size(Combos, 1);

    newBest = false;
    for n = 1:numCombos
      %Piece together which indices to use:
      ind = [bestInd, Combos(n,:)];
      
      %Now get the reduced sum of squares from the fit.
      
      [score_n, RSquared_n, coefs_n, yPredict_n] = ...
	  getLinearFit(x(:,ind), y_m);
      if score_n < rChiSquared(m)
	newBest = true;
	rChiSquared(m) = score_n;
	RSquared(m) = RSquared_n;
	coefs{m} = coefs_n;
	yPredictCell{m} = yPredict_n;
	fitIndices{m} = Combos(n,:);
      end
    end
    if newBest
      numAdd = 1;
      bestInd = fitIndices{m};
      if length(fitIndices{m}) == numX
	Done = true;
      end
    elseif(numAdd == options.maxLookAhead || ...
	   length(bestInd) + numAdd == numX || ...
	   length(bestInd) + numAdd == options.maxCombos)
      Done = true;
    else
      numAdd = numAdd + 1;
    end
  end
  
  if options.getPVals
    [YJack, pVals(m)] = JackknifeFit(x(:,bestInd), y_m);
  end
end

varargout = {fitIndices, yCell, yPredictCell, coefs};
if options.getPVals
  varargout = {varargout{:}, pVals};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rChiSquared, RSquared, varargout] = getLinearFit(x,y)
% [rChiSquared, coefs, yPredict] = getLinearFit(x, y)
numData = size(x, 1);
numXTypes = size(x, 2);
% numYTypes = size(y, 2);  %this equals 1

goodInd = find(sum(~isfinite(x), 2) == 0 & isfinite(y));
numGood = length(goodInd);
degreesOfFreedom = numGood - numXTypes - 1;
if(degreesOfFreedom < 1)
  rChiSquared = Inf;
  RSquared = 0;
  if nargout == 0
    varargout = {};
  elseif nargout == 1
    varargout = {repmat(NaN, numXTypes, 1)};
  else
    varargout = {repmat(NaN, numXTypes, 1), ...
		 repmat(NaN, numData, 1)};
  end
  return
end

yPredict = repmat(NaN, numData, 1);
x = [x, ones(numData, 1)];

coefs = pinv(x(goodInd,:)) * y(goodInd);
yPredict(goodInd) = x(goodInd,:) * coefs;

errY = yPredict(goodInd) - y(goodInd);
chiSquared = cov(errY) / cov(y(goodInd));
RSquared = 1 - chiSquared;
if RSquared < -1e-12 || chiSquared < -1e-12
  fprintf(2, 'Weird fit results.\n')
  keyboard
end
rChiSquared = chiSquared / degreesOfFreedom;

if nargout == 0
  varargout = {};
elseif nargout == 1
  varargout = {coefs};
else
  varargout = {coefs, yPredict};
end
return