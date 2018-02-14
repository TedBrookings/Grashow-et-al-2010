function [MapQuantities, MapQuantitiesLabels] = ...
    GetMapQuantities(Experiments, NetType)

MapQuantitiesLabels = {...
   'Frac HalfCenter', ...
   'Mean g_syn', 'Mean g_h', 'Mean f', ...
   'Max g_syn', 'Max g_h', 'Max f', ...
   'Min g_syn', 'Min g_h', 'Min f', ...
   'Span g_syn', 'Span g_h', 'Span f', ...
   'Std g_syn', 'Std g_h', 'Std f', ...
   'Mean AutoCorr', 'Mean SpikesPerBurst', 'Mean BurstSpikeFreq', ...
   'Mean DutyCycle', 'Mean SlowWaveAmp', ...
   'Max AutoCorr', 'Max SpikesPerBurst', 'Max BurstSpikeFreq', ...
   'Max DutyCycle', 'Max SlowWaveAmp', ...
   'Min AutoCorr', 'Min SpikesPerBurst', 'Min BurstSpikeFreq', ...
   'Min DutyCycle', 'Min SlowWaveAmp', ...
   'Span AutoCorr', 'Span SpikesPerBurst', 'Span BurstSpikeFreq', ...
   'Span DutyCycle', 'Span SlowWaveAmp', ...
   'Std AutoCorr', 'Std SpikesPerBurst', 'Std BurstSpikeFreq', ...
   'Std DutyCycle', 'Std SlowWaveAmp' ...
  };

              % Keep track of these:
IDList = {};              % unique experiment IDs
NumNetworks = [];         % number of networks in each experiment
StructList = [];          % some structure data
FList = {};               % frequencies
hList = {};               % g_h
synList = {};             % g_syn
CorrList = {};            % slow wave autocorrelation
SpikesPerBurstList = {};  % spikes per burst
BurstSpikeFreqList = {};      % ISI for spikes within a burst
DutyCycleList = {};       % Burst duty cycle
SlowWaveAmpList = {};     % Slow wave amplitude
NumHalfCenter = [];       % number of half centers

if(NetType == round(NetType))
  Flag = 'Category';
else
  Flag = 'AutoCorr';
  NetType = [NetType, Inf];
end

CellTypeList = GetCellTypeList(Experiments);
for n=1:length(CellTypeList)
  CellType = CellTypeList{n};
  
  %Get the networks of type CellType, with network properties NetType
  NetList = GetSpecifiedAnalysis(Experiments, 'CellType', CellType, ...
				 Flag, NetType);
  NumCells = length(unique({NetList.ID}));
  if(NumCells < 4)
    continue
  end
  %Get all networks of type CellType
  AllNetList = GetSpecifiedAnalysis(Experiments, 'CellType', CellType);
  
  for m = 1:length(NetList)
    %We are looping through all bursters of a given cell type ('GM', etc)

    %Get a few basic properties of this network
    NetStruct = NetList(m);
    syn_m = NetStruct.g_syn;  %g_syn
    h_m = NetStruct.g_h;      %g_h
    ID_m = NetStruct.ID;      %A unique identifying name for the experiment it's in
    Freq_m = NetStruct.CellReal.Burst.Freq;
    Corr_m = NetStruct.CellReal.SlowWave.Corr;
    SpikesPerBurst_m = NetStruct.CellReal.Burst.SpikesPerBurst.Mean;
    BurstSpikeFreq_m = NetStruct.CellReal.Burst.SpikeFreq;
    DutyCycle_m = NetStruct.CellReal.Burst.DutyCycle;
    SlowWaveAmp_m = mean(NetStruct.CellReal.SlowWave.Amplitudes);
    
    [MatchID, IDInd] = ismember(ID_m, IDList);
    if(~MatchID)  %if this is the first time we've seen this ID, make a new
                  %  "row' for this experiment
      IDList = {IDList{:}, ID_m};
      IDInd = length(IDList);
      Temp.CellType = CellType;
      Temp.ID = ID_m;
      Temp.Mean = 0;
      StructList = [StructList, Temp];
      
      NumNetworks = [NumNetworks, sum(strcmp({AllNetList.ID}, ID_m))];
      NumHalfCenter = [NumHalfCenter, 0];
      FList = {FList{:}, []};
      hList = {hList{:}, []};
      synList = {synList{:}, []};
      CorrList = {CorrList{:}, []};
      SpikesPerBurstList = {SpikesPerBurstList{:}, []};
      BurstSpikeFreqList = {BurstSpikeFreqList{:}, []};
      DutyCycleList = {DutyCycleList{:}, []};
      SlowWaveAmpList = {SlowWaveAmpList{:}, []};
    end
    
    %Now we know what "row" the experiment belongs in, so keep 
    %  "column" of map properties we're interested in
    NumHalfCenter(IDInd) = NumHalfCenter(IDInd) + 1;
    FList{IDInd} = [FList{IDInd}, Freq_m];
    hList{IDInd} = [hList{IDInd}, h_m];
    synList{IDInd} = [synList{IDInd}, syn_m];
    CorrList{IDInd} = [CorrList{IDInd}, Corr_m];
    SpikesPerBurstList{IDInd} = [SpikesPerBurstList{IDInd}, SpikesPerBurst_m];
    BurstSpikeFreqList{IDInd} = [BurstSpikeFreqList{IDInd}, BurstSpikeFreq_m];
    DutyCycleList{IDInd} = [DutyCycleList{IDInd}, DutyCycle_m];
    SlowWaveAmpList{IDInd} = [SlowWaveAmpList{IDInd}, SlowWaveAmp_m];
  end
end

%Calculate final map scalars
ListSize = size(NumHalfCenter);
FracHalfCenter = zeros(ListSize);

Mean_f = zeros(ListSize);
Mean_syn = zeros(ListSize);
Mean_h = zeros(ListSize);
Mean_AutoCorr = zeros(ListSize);
Mean_SpikesPerBurst = zeros(ListSize);
Mean_BurstSpikeFreq = zeros(ListSize);
Mean_DutyCycle = zeros(ListSize);
Mean_SlowWaveAmp = zeros(ListSize);

Max_f = zeros(ListSize);
Max_syn = zeros(ListSize);
Max_h = zeros(ListSize);
Max_AutoCorr = zeros(ListSize);
Max_SpikesPerBurst = zeros(ListSize);
Max_BurstSpikeFreq = zeros(ListSize);
Max_DutyCycle = zeros(ListSize);
Max_SlowWaveAmp = zeros(ListSize);

Min_f = zeros(ListSize);
Min_syn = zeros(ListSize);
Min_h = zeros(ListSize);
Min_AutoCorr = zeros(ListSize);
Min_SpikesPerBurst = zeros(ListSize);
Min_BurstSpikeFreq = zeros(ListSize);
Min_DutyCycle = zeros(ListSize);
Min_SlowWaveAmp = zeros(ListSize);

Std_f = zeros(ListSize);
Std_syn = zeros(ListSize);
Std_h = zeros(ListSize);
Std_AutoCorr = zeros(ListSize);
Std_SpikesPerBurst = zeros(ListSize);
Std_BurstSpikeFreq = zeros(ListSize);
Std_DutyCycle = zeros(ListSize);
Std_SlowWaveAmp = zeros(ListSize);

AutoCorr = zeros(ListSize);
SpikesPerBurst = zeros(ListSize);
BurstSpikeFreq = zeros(ListSize);
DutyCycle = zeros(ListSize);
SlowWaveAmp = zeros(ListSize);

for n = 1:length(FList)
  FracHalfCenter(n) = NumHalfCenter(n) / NumNetworks(n);

  [Mean_f(n), Max_f(n), Min_f(n), Span_f(n), Std_f(n)] = ExtractNums(FList{n});
  [Mean_syn(n), Max_syn(n), Min_syn(n), Span_syn(n), Std_syn(n)] = ...
      ExtractNums(synList{n});
  [Mean_h(n), Max_h(n), Min_h(n), Span_h(n), Std_h(n)] = ExtractNums(hList{n});

  [Mean_AutoCorr(n), Max_AutoCorr(n), Min_AutoCorr(n), Span_AutoCorr(n), ...
      Std_AutoCorr(n)] = ExtractNums(CorrList{n});
  [Mean_SpikesPerBurst(n), Max_SpikesPerBurst(n), Min_SpikesPerBurst(n), ...
      Span_SpikesPerBurst(n), Std_SpikesPerBurst(n)] ...
      = ExtractNums(SpikesPerBurstList{n});
  [Mean_BurstSpikeFreq(n), Max_BurstSpikeFreq(n), Min_BurstSpikeFreq(n), ...
      Span_BurstSpikeFreq(n), Std_BurstSpikeFreq(n)] = ...
      ExtractNums(BurstSpikeFreqList{n});
  [Mean_DutyCycle(n), Max_DutyCycle(n), Min_DutyCycle(n), ...
      Span_DutyCycle(n), Std_DutyCycle(n)] = ExtractNums(DutyCycleList{n});
  [Mean_SlowWaveAmp(n), Max_SlowWaveAmp(n), Min_SlowWaveAmp(n), ...
      Span_SlowWaveAmp(n), Std_SlowWaveAmp(n)] ...
      = ExtractNums(SlowWaveAmpList{n});
end

MapQuantitiesLabels = {...
   'Frac HalfCenter', ...
   'Mean g_syn', 'Mean g_h', 'Mean f', ...
   'Max g_syn', 'Max g_h', 'Max f', ...
   'Min g_syn', 'Min g_h', 'Min f', ...
   'Span g_syn', 'Span g_h', 'Span f', ...
   'Std g_syn', 'Std g_h', 'Std f', ...
   'Mean AutoCorr', 'Mean SpikesPerBurst', 'Mean BurstSpikeFreq', ...
   'Mean DutyCycle', 'Mean SlowWaveAmp', ...
   'Max AutoCorr', 'Max SpikesPerBurst', 'Max BurstSpikeFreq', ...
   'Max DutyCycle', 'Max SlowWaveAmp', ...
   'Min AutoCorr', 'Min SpikesPerBurst', 'Min BurstSpikeFreq', ...
   'Min DutyCycle', 'Min SlowWaveAmp', ...
   'Span AutoCorr', 'Span SpikesPerBurst', 'Span BurstSpikeFreq', ...
   'Span DutyCycle', 'Span SlowWaveAmp', ...
   'Std AutoCorr', 'Std SpikesPerBurst', 'Std BurstSpikeFreq', ...
   'Std DutyCycle', 'Std SlowWaveAmp' ...
  };

MapQuantities = ...
  {FracHalfCenter, ...
   Mean_syn, Mean_h, Mean_f, ...
   Max_syn, Max_h, Max_f, ...
   Min_syn, Min_h, Min_f, ...
   Span_syn, Span_h, Span_f, ...
   Std_syn, Std_h, Std_f, ...
   Mean_AutoCorr, Mean_SpikesPerBurst, Mean_BurstSpikeFreq, Mean_DutyCycle, ...
   Mean_SlowWaveAmp, ...
   Max_AutoCorr, Max_SpikesPerBurst, Max_BurstSpikeFreq, Max_DutyCycle, ...
   Max_SlowWaveAmp, ...
   Min_AutoCorr, Min_SpikesPerBurst, Min_BurstSpikeFreq, Min_DutyCycle, ...
   Min_SlowWaveAmp, ...
   Span_AutoCorr, Span_SpikesPerBurst, Span_BurstSpikeFreq, Span_DutyCycle, ...
   Span_SlowWaveAmp, ...
   Std_AutoCorr, Std_SpikesPerBurst, Std_BurstSpikeFreq, Std_DutyCycle, ...
   Std_SlowWaveAmp ...
  };


if(length(MapQuantities) ~= length(MapQuantitiesLabels))
  error('MapQuantities and Labels do not correspond.')
end
MapQuantities = AddStructInfo(MapQuantities, StructList);
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
function StripStr = StripNums(CellTypeStr)
StripStr = regexp(CellTypeStr, '[a-zA-Z]*', 'match');
if(length(StripStr) == 0)
  ErrStr = fprintf('Cell type %s is invalid:  %s', CellTypeStr, ...
		   'must contain letters');
  error(ErrStr)
else
  StripStr = StripStr{1};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MapQuantities = AddStructInfo(MapQuantities, StructList)
for n = 1:length(MapQuantities)
  TempList = MapQuantities{n};
  NewList = StructList;
  for m = 1:length(TempList)
    NewList(m).Mean = TempList(m);
  end
  MapQuantities{n} = NewList;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MeanVal, MaxVal, MinVal, Span, StdVal] = ExtractNums(List)
Ind = find(isfinite(List));
List = List(Ind);
if(length(List) > 2)
  MeanVal = mean(List);
  MaxVal = max(List);
  MinVal = min(List);
  Span = MaxVal - MinVal;
  StdVal = std(List);
else
  if(length(List) == 1)
    MeanVal = List;
  else
    MeanVal = NaN;
  end
  MaxVal = NaN;
  MinVal = NaN;
  Span = NaN;
  StdVal = NaN;
end
return