function DrawMeanChange(Organized)

ExpList = GetExperimentList(Organized);
Conditions = Organized.Conditions;
for n = 1:length(Conditions)
  Condition = Conditions{n};
  if(strcmp(Condition, 'ptx'))
    continue;
  end

  NoZero = 0;  %don't exclude zero-frequency from consideration
  for SpikesPerBurst = 0:2
    %if SpikesPerBurst = 0, analyze everything
    %if SpikesPerBurst = 1, analyze only networks with 1 spike per
    %                       burst in either condition
    %if SpikesPerBurst > 1, analyze only networks with spikes per
    %                       burst > 1 in both conditions
    [FreqControl, FreqMod] = GetAvgFreqs(Organized, ExpList, NoZero, ...
					 SpikesPerBurst, Condition);
    PlotMeanChangeFigure(FreqControl, FreqMod, NoZero, SpikesPerBurst, ...
			 Condition);
  end
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
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FreqControl, FreqMod] = GetAvgFreqs(Organized, ExpList, ...
					      ExcludeZeros, ...
					      SpikesPerBurst, Condition)
FreqControl = zeros(size(ExpList));
FreqMod = zeros(size(ExpList));
NumExp = length(ExpList);
Type_Mod = sprintf('Type_%s', Condition);
SpikesPerBurst_ptx = GetSpikesPerBurst(Organized.Analysis_ptx);
Analysis_mod = Organized.(sprintf('Analysis_%s', Condition));
SpikesPerBurst_mod = GetSpikesPerBurst(Analysis_mod);

for n = 1:NumExp
  ExpNum = ExpList(n);
  Ind = find(Organized.ExpNum == ExpNum);
  
  FListControl = Organized.Type_ptx(Ind);
  FListMod = Organized.(Type_Mod)(Ind);
  SPB_ptx = SpikesPerBurst_ptx(Ind);
  SPB_mod = SpikesPerBurst_mod(Ind);
  
  if(ExcludeZeros)
    ControlInd = find(FListControl ~= 0 & ...
                    isfinite(FListControl) & ...
                    isfinite(FListMod));
  else
    ControlInd = find(isfinite(FListControl) & ...
                    isfinite(FListMod));
  end
  switch(SpikesPerBurst)
   case 0,
    SpikeString = '';
   case 1,
    SpikeString = ' and One Spike Per Burst';
    Ind2 = find(SPB_ptx == 1 | SPB_mod == 1);
    ControlInd = intersect(ControlInd, Ind2);
   otherwise,
    SpikeString = ' and >1 Spike Per Burst';
    Ind2 = find(SPB_ptx > 1 & SPB_mod > 1);
    ControlInd = intersect(ControlInd, Ind2);
  end
  ModInd = ControlInd;
  
  FreqControl(n) = sum(FListControl(ControlInd)) / length(ControlInd);
  FreqMod(n) = sum(FListMod(ModInd)) / length(ModInd);
end

[TTestVal, PVal] = ttest(FreqControl, FreqMod, .05, 'both');
if(ExcludeZeros)
    ExcludeString = 'zeros excluded';
else
    ExcludeString = 'zeros not excluded';
end
TestString = sprintf('t-test with %s%s for %s:  p-value = %g', ...
    ExcludeString, SpikeString, Condition, PVal);
disp(TestString)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotMeanChangeFigure(FreqControl, FreqMod, ExcludeZeros, ...
			      SpikesPerBurst, Condition)
if(ExcludeZeros)
    ExcludeString = 'Zeros Excluded';
else
    ExcludeString = 'Zeros Not Excluded';
end
switch(SpikesPerBurst)
 case 0,
  SpikeString = '';
 case 1,
  SpikeString = ' and One Spike Per Burst';
 otherwise,
  SpikeString = ' and >1 Spike Per Burst';
end
  


Title = sprintf('Mean Change for %s with %s%s', Condition, ...
		ExcludeString, SpikeString);
h = NamedFigure(Title);

% setting colors
my_colors = colormap(cool(length(FreqControl)));
set(h,'DefaultAxesColorOrder',my_colors);
set(h, 'WindowStyle', 'docked');

xMat = repmat([1 2], length(FreqControl), 1);
yMat = [FreqControl', FreqMod'];

plot(xMat', yMat', '-o');
ylabel('Frequency (Hz)');
axis([0 3 -0.1 1.2]);
title(Title);
Axes_h = get(h, 'CurrentAxes');
set(Axes_h, 'FontName', 'Ariel');
set(Axes_h, 'FontSize', [15]);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SpikesPerBurst = GetSpikesPerBurst(AnalysisList)
Len = length(AnalysisList);
SpikesPerBurst = zeros(Len, 2);
for n=1:Len
  Temp = AnalysisList(n);
  if(length(Temp.Cell0) == 0 || length(Temp.Cell1) == 0)
    continue
  end
  try
    SpikesPerBurst(n,1) = Temp.Cell0.Burst.NumSpikes.Mean;
    SpikesPerBurst(n,2) = Temp.Cell1.Burst.NumSpikes.Mean;
    if(SpikesPerBurst(n,1) == 0)
      SpikesPerBurst(n,1) = 1;
    end
    if(SpikesPerBurst(n,2) == 0)
      SpikesPerBurst(n,2) = 1;
    end
  catch
    Temp
    whos
    error('poo')
  end
end

SpikesPerBurst = max(SpikesPerBurst, [], 2);
return