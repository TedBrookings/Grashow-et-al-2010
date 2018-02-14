function XR = ZScore(X, Dim)
if(nargin < 2)
  Dim = 1;
  while(size(X, Dim) < 2)
    if(size(X, Dim) == 0)
      XR = X;
      return
    end
    Dim = Dim + 1;
    if(Dim > ndims(X))
      XR = X;
      return
    end
  end
end

NonFinite = (sum(~isfinite(X(:))) > 0);
if(NonFinite)
  XR = repmat(NaN, size(X));
  switch(Dim)
   case 1,
    for m = 1:size(X, 2);
      GoodInd = find(isfinite(X(:,m)));
      XR(GoodInd,m) = zscore(X(GoodInd,m));
    end
   case 2,
    for m = 1:size(X, 1);
      GoodInd = find(isfinite(X(m,:)));
      XR(m,GoodInd) = zscore(X(m,GoodInd));
    end    
   otherwise,
    error('RankScore only works on matrices of dimension <= 2')
  end    
else
  XR = zscore(X, Dim);
end
return