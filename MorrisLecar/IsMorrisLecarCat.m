function Test = IsMorrisLecarCat(Cat, Test)
if(nargin < 2)
  Test = 0;
end

CatList = {'5ht', 'oxo', 'da', 'dc'};
MLType = true;
for n=1:length(CatList)
  TestCat = CatList{n};
  %whos
  %Str = sprintf('%s, %s:  %g', Cat, TestCat, length(regexp(Cat, TestCat)));
  %disp(Str)
  if(length(regexp(Cat, TestCat)) > 0)
    MLType = false;
    break;
  end
end

if(MLType)
  if(Test == 0)
    Test = 1;
  elseif(Test == -1)
    error(['Mixture of Morris-Lecar experiments and others within' ...
	   ' one directory!'])
  end
else
  if(Test == 0)
    Test = -1;
  elseif(Test == 1)
    error(['Mixture of Morris-Lecar experiments and others within' ...
	   ' one directory!'])
  end
end
return