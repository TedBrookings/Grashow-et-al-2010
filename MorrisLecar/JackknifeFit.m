function [YPredict, varargout] = JackknifeFit(X, Y)
% [YPredict, PVal, RSquared] = JackknifeFit(X,Y)
% Performs jackknifed linear fit of Y to X:  Y = A * X + B
%  INPUTS:
%   -X: NData x NX-types matrix
%   -Y: NData x NY-types matrix
%  OUTPUTS:
%   -YPredict: NData x NY-types matrix of predicted Y values.
%               The nth YPredict value is predicted from a linear
%               fit that *excludes* the nth X and Y values.
%   OPTIONAL:
%   -PVal:  NY-types array of p-values, calculated from Pearson
%               correlation value between Y and YPredict.
%   -RSquared: NY-types array of R^2 values, between Y and YPredict.

if(nargout > 3)
  error('Too many output arguments.  Run "help JackknifeFit"')
end
NData = size(X, 1);
if(size(Y,1) ~= NData)
  error('X and Y must have the same number of rows')
end
NYTypes = size(Y, 2);
YPredict = repmat(NaN, NData, NYTypes);
X = [X, ones(NData,1)];

GoodX = find(sum(~isfinite(X),2) == 0);
NGoodX = length(GoodX);

for n = 1:NGoodX
  Ind = GoodX([1:(n-1),(n+1):NGoodX]);
  m = GoodX(n);
  for k = 1:NYTypes
    Ind_k = Ind(find(isfinite(Y(Ind,k))));
    Coef = pinv(X(Ind_k,:)) * Y(Ind_k,k);
    YPredict(m,k) = X(m,:) * Coef;
  end
end

switch(nargout)
 case 1,
  varargout = {};
 case 2,
  PVal = zeros(NYTypes,1);
  for n = 1:NYTypes
    [Coef, PVal(n)] = Pearson(Y(:,n), YPredict(:,n));
  end
  varargout = {PVal};
 case 3,
  PVal = zeros(NYTypes,1);
  RSquared = zeros(NYTypes,1);
  for n = 1:NYTypes
    [Coef, PVal(n)] = Pearson(Y(:,n), YPredict(:,n));
    FiniteInd = find(isfinite(Y(:,n)) & isfinite(YPredict(:,n)));
    
    YDev = Y(FiniteInd,n);
    YErr = YPredict(FiniteInd,n) - YDev;
    %YDev = YDev - mean(YDev);
    %RSquared(n) = 1 - YErr' * YErr / (YDev' * YDev);
    RSquared(n) = 1 - cov(YErr) / cov(YDev);
    if(isnan(RSquared(n)))
      fprintf(2, 'NaN in JackknifeFit.m\n')
      keyboard
    end
  end
  varargout = {PVal, RSquared};
end

return