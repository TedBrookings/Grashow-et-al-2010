function XR = RankScore(X, Dim)
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

NonFinite = (sum(isnan(X(:))) > 0);
if(NonFinite)
  XR = repmat(NaN, size(X));
  switch(Dim)
   case 1,
    for m = 1:size(X, 2);
      GoodInd = find(~isnan(X(:,m)));
      [XSort, Ind] = sort(X(GoodInd,m));
      IndArr = linspace(0, 1, length(Ind));
      XR(GoodInd(Ind),m) = IndArr;
    end
   case 2,
    for m = 1:size(X, 1);
      GoodInd = find(~isnan(X(m,:)));
      [XSort, Ind] = sort(X(m,GoodInd));
      IndArr = linspace(0, 1, length(Ind));
      XR(m,GoodInd(Ind)) = IndArr;
    end    
   otherwise,
    error('RankScore only works on matrices of dimension <= 2')
  end    
else
  [XSort, Ind] = sort(X, Dim);
  Len = size(X, Dim);
  IndArr = linspace(0, 1, Len);
  switch(Dim)
   case 1,
    for m = 1:size(X, 2)
      XR(Ind(:,m), m) = IndArr;
    end
   case 2,
    for m = 1:size(X, 1);
      XR(m, Ind(m,:)) = IndArr;
    end
   otherwise,
    error('RankScore only works on matrices of dimension <= 2')
  end
end
return