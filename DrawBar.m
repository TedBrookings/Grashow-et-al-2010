function DrawBar(Organized)
MakePairwisePlots = true;
MakeOverallPlots = true;

Conditions = Organized.Conditions;
NumConditions = length(Conditions);
BurstNums = [];
n_ptx = -1;

for n = 1:NumConditions
  Condition = Conditions{n};
  BurstNums = [BurstNums, GetBurstNums(Organized, Conditions{n})];

  if(strcmp(Conditions{n}, 'ptx'))
    n_ptx = n;
  end
end

if(MakePairwisePlots)
  for n=1:NumConditions
    if(n == n_ptx)
      continue;
    end
    
    MakeBurstBarPlot(Conditions, BurstNums, n);
  end
end

if(MakeOverallPlots)
  CondList = [1:(n_ptx-1), (n_ptx+1):NumConditions];
  MakeBurstBarPlot(Conditions, BurstNums, CondList);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BurstNums = GetBurstNums(Organized, Condition)
TypeStr = sprintf('Type_%s', Condition);
Values = Organized.(TypeStr);
GoodInd = find(isfinite(Values));
Values = Values(GoodInd);
ExpNum = Organized.ExpNum(GoodInd);
PtxValues = Organized.Type_ptx(GoodInd);

BurstNums.ExpList = unique(ExpNum);
NumExp = length(BurstNums.ExpList);
BurstNums.NumTotal = zeros(1, NumExp);
BurstNums.NumBurst = zeros(1, NumExp);
BurstNums.NumPtx = zeros(1, NumExp);
for n=1:NumExp
  ExpInd = find(ExpNum == BurstNums.ExpList(n));
  BurstNums.NumTotal(n) = length(ExpInd);
  BurstNums.NumBurst(n) = length(find(Values(ExpInd) > 0));
  BurstNums.NumPtx(n) = length(find(PtxValues(ExpInd) > 0));
end

if(strcmp(Condition, 'ptx'))
  return
end
[Sig, pVal] = ttest(BurstNums.NumBurst - BurstNums.NumPtx);
disp(sprintf(['T-test for number of bursters for %s vs ptx', ...
	      ' has p=%g'], Condition, pVal))
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeBurstBarPlot(Conditions, BurstNums, IndList)
PVal = erfc(1/sqrt(2));  %1 sigma
Num = length(IndList);
BurstFrac = zeros(Num, 1);
BurstStd = zeros(Num, 1);
PtxFrac = zeros(Num, 1);
PtxStd = zeros(Num, 1);
y = zeros(2 * Num, 1);
yConf = zeros(2 * Num, 2);
Labels = {};
for n=1:Num
  Total = BurstNums(n).NumTotal;
  Burst = BurstNums(n).NumBurst ./ Total;
  Ptx = BurstNums(n).NumPtx ./ Total;
  
  BurstFrac(n) = mean(Burst);
  BurstStd(n) = std(Burst);
  PtxFrac(n) = mean(Ptx);
  PtxStd(n) = std(Ptx);
  
  y(2*n - 1) = PtxFrac(n);
  yConf(2*n - 1, 1) = PtxFrac(n) - PtxStd(n);
  yConf(2*n - 1, 2) = PtxFrac(n) + PtxStd(n);
  y(2*n) = BurstFrac(n);
  yConf(2*n,1) = BurstFrac(n) - BurstStd(n);
  yConf(2*n,2) = BurstFrac(n) + BurstStd(n);
  Labels = {Labels{:}, [Conditions{IndList(n)}, ' control'], ...
	    Conditions{IndList(n)}};
end

if(Num == 1)
  TitleStr = sprintf('Overall Bursters for ptx and %s', ...
		     Conditions{IndList});
  yStr = sprintf('Overall Percent of Bursters for ptx and %s', ...
		 Conditions{IndList});
else
  TitleStr = 'Overall Bursters';
  yStr = 'Overall Percent of Bursters';
end

h = NamedFigure(TitleStr);
set(h, 'WindowStyle', 'docked');
clf

x = 1:(2*Num);
bar(100 * y);
hold on
for n = x
  plot([n, n], yConf(n,:), 'r-');
end
hold off

Axes_h = get(h, 'CurrentAxes');
set(Axes_h, 'XTick', x);
set(Axes_h, 'XTickLabel', Labels);
xlabel('Condition');
ylabel(yStr);
title(TitleStr);
set(Axes_h, 'FontName', 'Ariel');
set(Axes_h, 'FontSize', [15]);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Organized = GetDecrease(Organized, ChangeThreshold, Condition)
CondStr = sprintf('Type_%s', Condition);
Organized.Change = Organized.(CondStr) - Organized.Type_ptx;
Ind = find(Organized.Change < -ChangeThreshold);
Names = fieldnames(Organized);
for n = 1:length(Names)
  FieldName = Names{n};
  if(strcmp(FieldName, 'Conditions'))
      continue;
  end
  Arr = getfield(Organized, FieldName);
  Arr = Arr(Ind);
  Organized = setfield(Organized, FieldName, Arr);
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExpList = GetExperimentList(Organized)
Nums = Organized.ExpNum;
ExpList = Nums(1);
for n = 2:length(Nums)
  if(length(find(ExpList == Nums(n))) == 0)
    ExpList = [ExpList, Nums(n)];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OldMakeExperimentPlot(Organized, ExpList, Condition)

ExpNum = Organized.ExpNum;
Total = zeros(size(ExpList));
for n = 1:length(ExpList)
    Total(n) = length(find(ExpNum == ExpList(n)));
end

TitleStr = sprintf('Networks Decreasing Under %s', Condition);

h = NamedFigure(TitleStr);
set(h, 'WindowStyle', 'docked');
clf;

bar(Total);
set(gca, 'XTick', 1:length(Total));
set(gca, 'XTickLabel', ExpList);
xlabel('Experiment Number');
ylabel(['Number of ', TitleStr]);
title(TitleStr);
return
