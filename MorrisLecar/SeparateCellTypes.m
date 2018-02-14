%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SeparateCellTypes(IPs, mapProps, mapList, options)

progressiveIP = true;
progressiveMapProps = true;

cellTypes = IPs.cellType;

if(progressiveIP)
  bestIP = progressiveSeparateTypes(IPs.mat, IPs.labels, ...
				    cellTypes, 'IP');
  IPs.mat = IPs.mat(:, bestIP);
else
  separateTypes(IPs.mat, cellTypes, 'IP');
end

if(progressiveMapProps)
  bestMapProps = progressiveSeparateTypes(mapProps.mat, ...
					  mapProps.labels, cellTypes, ...
					  'Map Prop');
  mapProps.mat = mapProps.mat(:, bestMapProps);
else
  separateTypes(mapProps.mat, cellTypes, 'Map Prop');
end

separateTypes(mapList, cellTypes, 'Map Pos', @mapDistance);

separateTypes([IPs.mat, mapProps.mat], cellTypes, 'IP + Map Prop');
separateTypes({IPs.mat, mapList}, cellTypes, 'IP + Map Pos', ...
	      {'Euclidean', @mapDistance});
separateTypes({mapProps.mat, mapList}, cellTypes, ...
	      'Map Prop + Map Pos', {'Euclidean', @mapDistance});

separateTypes({[IPs.mat, mapProps.mat], mapList}, cellTypes, ...
	      'IP + Map Prop + Map Pos', {'Euclidean', @mapDistance});

separateTypes(randn(size(mapProps.mat)), cellTypes, 'Random Nums');
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = separateTypes(dataMat, cellTypes, label, varargin)
%numClustTrials = 1000;
numClustTrials = Inf;  %does exhaustive search

%Remove cells that have all NaN properties
try
  okayCells = find(sum(isnan(dataMat), 2) < size(dataMat, 2));
  dataMat = dataMat(okayCells,:);
  cellTypes = cellTypes(okayCells);
catch dataNotDoubles
  ;
end

[typeList, lastInds, typeInds] = unique(cellTypes);
numTypes = length(typeList);
numCells = length(cellTypes);
clustInds = kmedoids(dataMat, numTypes, numClustTrials, varargin{:});
typeInds = typeInds';

%Fill up a matrix of what clusters each type was assigned to.
indMat = zeros(numTypes, numTypes);
for n = 1:numTypes
  ind_n = find(typeInds == n);
  clust_n = clustInds(ind_n);
  hist_clust_n = histc(clust_n, 0.5 + (0:numTypes));
  indMat(n,:) = hist_clust_n(1:numTypes);
end

%Swap around as necessary, to assign the right labels to each
%cluster
for n = 1:numTypes
  [maxVal, maxCol] = max(indMat(n,:));
  if(indMat(n,n) == maxVal)
    continue
  end
  if(indMat(n,n) + indMat(maxCol,maxCol) > ...
     indMat(n,maxCol) + indMat(maxCol,n))
    continue
  end
  %swap columns
  temp = indMat(:,n);
  indMat(:,n) = indMat(:,maxCol);
  indMat(:,maxCol) = temp;
end

fracRight = trace(indMat) / numCells;

%save results in compact structure:
clustResults.label = label;
clustResults.fracRight = fracRight;
clustResults.indMat = indMat;
clustResults.typeList = typeList;

%display results
if(nargout == 0)
  displayClusterResults(clustResults);
elseif(size(dataMat, 2) == 1)
  fprintf('Clustering by %s:\n', label)
  fprintf('  ->Fraction correct:  %g\n', fracRight)
end

if(nargout == 0)
  varargout = {};
else
  varargout = {clustResults};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayClusterResults(clustResults)
fprintf('Clustering by %s:\n', clustResults.label)
numTypes = length(clustResults.typeList);
for n = 1:numTypes
  fprintf('%s:', clustResults.typeList{n})
  for m = 1:numTypes
    fprintf(' %g', clustResults.indMat(n,m))
  end
  fprintf('\n')
end
fprintf('  ->Fraction correct:  %g\n', clustResults.fracRight)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bestResults, bestInd] = bestSeparateTypes(X, XLabels, ...
						  cellTypes, label, ...
						  numUse, varargin)
[typeList, lastInds, typeInds] = unique(cellTypes);
numTypes = length(typeList);
numXTypes = size(X, 2);

ind = nchoosek(1:numXTypes, numUse);
numTrials = size(ind, 1);
bestResults.fracRight = -Inf;
for n = 1:numTrials
  ind_n = ind(n,:);
  if(numUse == 1)
    label_n = [label, ' (', XLabels{ind_n}, ')'];
  else
    label_n = label;
  end
  clustResults = separateTypes(X(:,ind_n), cellTypes, label_n, ...
			       varargin{:});
  if(clustResults.fracRight > bestResults.fracRight)
    bestResults = clustResults;
    bestInd = ind_n;
  end
end

bestLabel = XLabels{bestInd(1)};
for k = 2:length(bestInd)
  bestLabel = [bestLabel, ', ', XLabels{bestInd(k)}];
end

fprintf('Clustering by %s: (%s)\n', label, bestLabel)
fprintf('  ->Fraction correct:  %g\n', bestResults.fracRight)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bestInd = progressiveSeparateTypes(dataMat, DataLabels, ...
					    cellTypes, label, varargin)
bestResults.fracRight = -Inf;
numTypes = size(dataMat, 2);
fprintf('Finding best separation for %s with %u types.\n', ...
	label, numTypes)
for numUse = 1:numTypes
  if(numUse == 1)
    label_n = sprintf('Best %u %s', numUse, label);
  else
    label_n = sprintf('Best %u %ss', numUse, label);
  end
  [clustResults, ind_n] = bestSeparateTypes(dataMat, DataLabels, ...
					    cellTypes, label_n, ...
					    numUse, varargin{:});
  if(clustResults.fracRight > bestResults.fracRight)
    bestResults = clustResults;
    bestInd = ind_n;
  end
end
label_n = sprintf('All %ss', label);
clustResults = separateTypes(dataMat, cellTypes, label_n, varargin{:});
if(clustResults.fracRight > bestResults.fracRight)
  bestResults = clustResults;
  bestInd = 1:numTypes;
end

fprintf('Best clustering:')
for n = 1:length(bestInd)
  fprintf(' %s', DataLabels{bestInd(n)})
end
fprintf('\n')
displayClusterResults(bestResults);
return