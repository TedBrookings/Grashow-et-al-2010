function [clustInds, varargout] = kmedoids(dataMat, k, numTries, ...
					   varargin)
% [clustInds, centroids] = kmedoids(dataMat, k, numTries, distType)
% Implements k-medoids++ clustering algorithm (k-medoids with a
%   particular initialization choice).
%  INPUTS:
%    dataMat: NumData x DataDim matrix of points
%             optionally, dataMat can be a cell array, with each
%               element having the same number of rows.  (This
%               option allows specification of different distance
%               functions, as well as the use of classes.)  When
%               computing distances the distance from each cell
%               is added to the total distance.
%    k:  Number of clusters to find
%   (OPTIONAL)
%    numTries: Take the best cluster scheme after this many tries.
%              (defaults to 1).  The best scheme is determined by
%              minimizing the average minimum distance from a
%              point to its centroid.
%    distType: Specify the function to determine distance between
%              two points.  Can pass a function handle, or specify
%              the following strings:
%               'Euclidean'
%              (defaults to 'Euclidean')
%              If dataMat is a cell array, distType should either
%              be empty, or a cell array of the same length.
%  OUTPUTS:
%    clustInds: NumData array of indices indicated which cluster
%               each point is assigned to.
%   (OPTIONAL)
%    centroids: k x DataDim matrix of cluster centroids

if(nargin < 3)
  numTries = 1;
end

if numTries < 1
  help kmedoids
  error('Invalid numTries:  %g\n', numTries)
end

[hClusterDist, dataMat] = getHClusterDist(dataMat, varargin{:});

if numTries >= size(dataMat, 1)
  numTries = size(dataMat, 1);
  exhaustiveStart = true;
else
  exhaustiveStart = false;
end

if numTries == 1
  [clustInds, meanDistance, centroids] = ...
      do_kmedoids(dataMat, k, hClusterDist);
else
  if exhaustiveStart
    ind = {1};
  else
    ind = {};
  end
  
  [clustInds, meanDistance] = do_kmedoids(dataMat, k, hClusterDist, ...
					  ind{:});
  for n = 2:numTries
    if exhaustiveStart
      ind = {n};
    end
    [clustInds_n, meanDistance_n, centroids_n] = ...
	do_kmedoids(dataMat, k, hClusterDist, ind{:});
    if meanDistance_n < meanDistance
      clustInds = clustInds_n;
      meanDistance = meanDistance_n;
      centroids = centroids_n;
    end
  end
end
if nargout == 2
  varargout = {centroids};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hClusterDist, dataMat] = getHClusterDist(dataMat, varargin)
funcNames = {'Euclidean', 'Taxi'};
defaultFunc = 'Euclidean';
funcHandles = {@euclideanDist, @taxiDist};

if iscell(dataMat)
  %there are multiple independent data sets, using %@multiDimCellDist
  numSets = length(dataMat);
  fHands = cell(1, numSets);
  if isempty(varargin)
    %use default
    distTypes = cell(1, numSets);
    for n = 1:numSets
      distTypes{n} = defaultFunc;
    end
  else
    %use user-specified
    distTypes = varargin{1};
  end
  
  for n = 1:numSets
    %get the handles
    fHands{n} = getFHandle(distTypes{n}, funcNames, funcHandles);
  end
  %wrap them
  hClusterDist = @(p1, p2) multiDimCellDist(p1, p2, fHands);
  
  %now reform dataMat to be useful in this form
  oldDataMat = dataMat;
  numRows = size(oldDataMat{1},1);
  dataMat = cell(numRows, numSets);
  for m = 1:numSets
    mat_m = oldDataMat{m};
    for n = 1:numRows
      dataMat{n,m} = mat_m(n,:);
    end
  end
else
  %just one data set
  if isempty(varargin)
    %use default
    distType = defaultFunc;
  else
    %use user-specified
    distType = varargin{1};
  end
  %get the handle
  hClusterDist = getFHandle(distType, funcNames, funcHandles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fHandle = getFHandle(distType, funcNames, funcHandles)
switch class(distType)
 case 'char',
  for n = 1:length(funcNames)
    if StringCheck(funcNames{n}, distType)
      fHandle = funcHandles{n};
      return
    end
  end
  error('Invalid distType:  %s\n', distType)
 case 'function_handle',
  fHandle = distType;
 otherwise,
  error('Invalid distType class:  %s\n', class(distType))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clustInds, varargout] = do_kmedoids(dataMat, k, hClusterDist, ...
					      varargin)

centroids = getFirstCentroids(dataMat, k, hClusterDist, varargin{:});

clustInds = zeros(size(dataMat, 1), 1);
[clustInds, meanDistance] = updateClustInds(dataMat, ...
					clustInds, centroids, ...
					hClusterDist);
if length(unique(clustInds)) ~= k
  fprintf(2, 'Some clusters not assigned at startup!\n')
  keyboard
end

converged = false;
numItr = 0;
maxItr = size(dataMat,1)^2;
while ~converged && numItr < maxItr
  centroids = updateCentroids(dataMat, centroids, clustInds, hClusterDist);

  oldClustInds = clustInds;
  [clustInds, meanDistance] = updateClustInds(dataMat, ...
					  clustInds, centroids, ...
					  hClusterDist);
  if length(unique(clustInds)) ~= k
    fprintf(2, 'Some clusters not assigned in step %d!\n', ...
	    numItr + 1)
    keyboard
  end
  
  numItr = numItr + 1;
  converged = (sum(clustInds ~= oldClustInds) == 0);
end
if nargout == 2
  varargout = {meanDistance};
elseif nargout == 3
   varargout = {meanDistance, centroids};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function centroids = getFirstCentroids(dataMat, k, hClusterDist, startInd)
numPoints = size(dataMat, 1);
clustDim = size(dataMat, 2);
centroidInds = zeros(k, 1);

%Choose all cluster means randomly by the following method:
%  -Choose the first point uniformly
%  -Choose each subsequent point from remaining points, with
%   probability proportional to distance from the nearest centroid

if nargin == 4
  chooseInd = startInd;
else
  chooseInd = 1 + floor(numPoints * rand(1,1));
end

n = 1;
centroidInds(n) = chooseInd;
probSelect = zeros(numPoints, 1);
centroids = repmat(dataMat(chooseInd,:), k, 1);

for n = 2:k
  %find the nth centroid, starting at n=2
  
  for ind = setdiff(1:numPoints, centroidInds(1:(n-1)))
    %for each non-centroids point, find the distance to the closest
    %  centroid
    minDistance = Inf;
    data_ind = dataMat(ind,:);
    for centroidInd = centroidInds(1:(n-1))'
      %loop over centroids, finding distance and comparing to min
      dist = hClusterDist(dataMat(centroidInd,:), data_ind);
	
      if dist < minDistance
	minDistance = dist;
      end
    end
    %the probability of selection is proportional to min distance
    probSelect(ind) = minDistance;
  end
  % make the probSelect cumulative
  probSelect = cumsum(probSelect);
  % choose proprortionally to probability
  chooseInd = find(probSelect >= rand(1,1) * probSelect(end), 1);
  
  centroidInds(n) = chooseInd;
  centroids(n,:) = dataMat(chooseInd,:);
  if(n < k)
    probSelect(:) = 0;
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clustInds, meanDistance] = updateClustInds(dataMat, ...
						 clustInds, centroids, ...
						 hClusterDist)
meanDistance = 0;
k = size(centroids, 1);
numPoints = size(dataMat, 1);

for n = 1:numPoints
  data_n = dataMat(n,:);
  minDist = Inf;
  for m = 1:k
    dist = hClusterDist(data_n, centroids(m,:));
    if dist < minDist
      minDist = dist;
      clustInds(n) = m;
    end
  end
  if ~ismember(clustInds(n), 1:k)
    fprintf(2, 'WTF clustInd(n) = %g\n', clustInds(n))
    keyboard
  end
  meanDistance = meanDistance + minDist;
end

meanDistance = meanDistance / numPoints;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function centroids = updateCentroids(dataMat, centroids, clustInds, ...
				     hClusterDist)
k = size(centroids, 1);

%loop through each cluster
for n = 1:k
  %get the data that are in this cluster
  ind = find(clustInds == n);
  num_n = length(ind);
  %try making each member of the cluster the centroid
  %  chose the centroid with the lowest within-cluster distance
  
  minDist = Inf;
  for m = 1:num_n
    dist_m = 0;
    data_m = dataMat(ind(m),:);
    for k = [1:(m-1), (m+1):num_n]
      data_k = dataMat(ind(k),:);
      dist_m = dist_m + hClusterDist(data_m, data_k);
    end
    if dist_m < minDist
      centroids(n,:) = data_m;
      minDist = dist_m;
    end
  end
end

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = euclideanDist(p1, p2)
finiteInd = find(isfinite(p1) & isfinite(p2));
diffP = p1(finiteInd) - p2(finiteInd);
dist = (diffP * diffP') / length(finiteInd);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = multiDimCellDist(p1, p2, hDist)
dist = hDist{1}(p1{1}, p2{1});
for n = 2:length(p1)
  dist = dist + hDist{n}(p1{n}, p2{n});
end
return