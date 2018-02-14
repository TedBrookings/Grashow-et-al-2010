function [synList, hList] = GetSensitivity(Experiments, NetType)

FreqType = 'Burst';
if(NetType == round(NetType))
  Flag = 'Category';
  if(NetType == 1)
    FreqType = 'Spike';
  end
else
  Flag = 'AutoCorr';
  NetType = [NetType, Inf];
end



CellTypeList = GetCellTypeList(Experiments);
synList = [];
hList = [];

for n=1:length(CellTypeList)
  CellType = CellTypeList{n};
  %Get the networks that are bursters of type CellType:
  NetList = GetSpecifiedAnalysis(Experiments, 'CellType', CellType, ...
				 Flag, NetType);
  ID = {NetList.ID};
  NumCells = length(unique(ID));
  if(NumCells < 4)
    continue
  end
  g_syn = cat(1, NetList.g_syn);
  g_h = cat(1, NetList.g_h);
  
  syn_Sensitivity.CellType = CellType;
  syn_Sensitivity.ID = [];
  syn_Sensitivity.Sens = [];
  syn_Sensitivity.g_syn = [];
  syn_Sensitivity.g_h = [];
  
  h_Sensitivity.CellType = CellType;
  h_Sensitivity.ID = [];
  h_Sensitivity.Sens = [];
  h_Sensitivity.g_syn = [];
  h_Sensitivity.g_h = [];

  %Doesn't work with sub-sub-structures:
  %Freq = cat(1, NetList.CellReal.Burst.Freq);
  NumNetworks = length(g_syn);
  for m = 1:NumNetworks
    g_syn_m = g_syn(m);
    g_h_m = g_h(m);
    Freq_m = NetList(m).CellReal.(FreqType).Freq;
    ID_m = NetList(m).ID;
    
    syn_Match = find(strcmp(ID, ID_m)' & g_h == g_h_m & g_syn == g_syn_m + 15);    

    h_Match = find(strcmp(ID, ID_m)' & g_syn == g_syn_m & g_h == g_h_m + 15);

    if(length(syn_Match) == 1)
      FreqMatch = NetList(syn_Match).CellReal.(FreqType).Freq;
      %Relative frequency change / Delta_g:
      Sens = 2.0 * (FreqMatch - Freq_m) / (FreqMatch + Freq_m) / 15;
      if(isnan(Sens))
        Sens = 0;
      end
      syn_Sensitivity.ID = ID_m;
      syn_Sensitivity.Sens = Sens;
      syn_Sensitivity.g_syn = 0.5 * (g_syn_m + g_syn(syn_Match));
      syn_Sensitivity.g_h = g_h(syn_Match);
      synList = [synList, syn_Sensitivity];
    elseif(length(syn_Match) > 1)
      error('Too many sensitivity matches.')
    end
    if(length(h_Match) == 1)
      FreqMatch = NetList(h_Match).CellReal.(FreqType).Freq;
      %Relative frequency change / Delta_g:
      Sens = 2.0 * (FreqMatch - Freq_m) / (FreqMatch + Freq_m) / 15;
      if(isnan(Sens))
        Sens = 0;
      end
      h_Sensitivity.ID = ID_m;
      h_Sensitivity.Sens = Sens;
      h_Sensitivity.g_syn = g_syn(h_Match);
      h_Sensitivity.g_h = 0.5 * (g_h_m + g_h(h_Match));
      hList = [hList, h_Sensitivity];
    elseif(length(h_Match) > 1)
      error('Too many sensitivity matches.')
    end
  end
end

%PlotSens(synList, hList);

synList = AverageOverCell(synList);
hList = AverageOverCell(hList);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CellTypeList = GetCellTypeList(Experiments)
CellTypeList = {};
for n = 1:length(Experiments)
  Temp = Experiments(n).ConditionList;
  for m = 1:length(Temp)
    Temp{m} = StripNums(Temp{m});
  end
  CellTypeList = {CellTypeList{:}, Temp{:}};
end
CellTypeList = unique(CellTypeList);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StripStr = StripNums(ConditionStr)
StripStr = regexp(ConditionStr, '[a-zA-Z]*', 'match');
if(length(StripStr) == 0)
  ErrStr = fprintf('Experiment type %s is invalid:  %s', ConditionStr, ...
		   'must contain letters');
  error(ErrStr)
else
  StripStr = StripStr{1};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotSens(synList, hList)
NetList = unique({synList.CellType, hList.CellType});

for n = 1:length(NetList)
  CellType = NetList{n};
  synInd = find(strcmp({synList.CellType}, CellType));
  hInd = find(strcmp({hList.CellType}, CellType));
  
  Title = sprintf('Sensitivity for %s', CellType);

  NumBins = round(length(synInd) / 5);
  if(NumBins < 10)
    NumBins = 10;
  end

  [synNum, synSens] = hist(cat(1,synList(synInd).Sens), NumBins);
  
  NumBins = round(length(hInd) / 5);
  if(NumBins < 10)
    NumBins = 10;
  end
  [hNum, hSens] = hist(cat(1,hList(hInd).Sens), NumBins);
  
  h = NamedFigure(Title);
  set(h, 'WindowStyle', 'docked');
  hold off
  bar(synSens, synNum, 'b')
  hold on
  bar(hSens, hNum, 'r')
  hold off
  xlabel('Relative sensitivity', 'FontName', 'Arial', 'FontSize', 18);
  ylabel('Number of networks', 'FontName', 'Arial', 'FontSize', 18);
  title(Title, 'FontName', 'Arial', 'FontSize', 18);
  Axes_h = get(h, 'CurrentAxes');
  set(Axes_h, 'FontName', 'Arial');
  set(Axes_h, 'FontSize', [30]);
  legend('g_s_y_n', 'g_h', 'Location', 'Best')
  xlim([-.15, .15])
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AvgSensList = AverageOverCell(SensList)
IDList = unique({SensList.ID});

AvgSensList = [];
for n = 1:length(IDList);
  Temp.ID = IDList{n};
  Ind = find(strcmp({SensList.ID}, Temp.ID));
  if(length(Ind) == 0)
    error('No cells with ID;  %s', Temp.ID)
  elseif(length(Ind) == 1)
    Temp.Mean = SensList(Ind).Sens;
    Temp.Std = NaN;
  else
    Sens = cat(2, SensList(Ind).Sens);
    Temp.Mean = mean(Sens);
    Temp.Std = std(Sens);
  end
  if(isnan(Temp.Mean) || length(Temp.Std) == 0)
    fprintf(2, 'Sensitivity for %s: Mean=%g, Std=%g\n', Temp.ID, Temp.Mean, ...
	                                                 Temp.Std)
    keyboard                                                        
  end
  AvgSensList = [AvgSensList, Temp];
end
return