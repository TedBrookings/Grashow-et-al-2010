function FI = GetFICurve(varargin)
% FI = GetFICurve(FileName, CellName, PlotVar, tCutoff)
%        ---OR---
% FI = GetFICurve(t, v, I, PlotVar, tCutoff)
% if no parameters are passed, a structure filled with NaNs is returned.
%  INPUTS:
%   -t is time in ms
%   -v is voltage in mV
%   -I is current in nA
%   -FileName is the name of an .abf file with FI data
%   -CellName is the name of the requested cell within the
%      file specified by FileName (e.g. "top" or "bot")
%    OPTIONAL:
%     -PlotVar set to true/1 to plot the result
%     (defaults to false/0)
%     -tCutoff sets time (in s) during current injection to look at spikes
%        if tCutoff > 0, looks at first tCutoff seconds of data
%        if tCutoff < 0, looks at last abs(tCutoff) seconds of data
%        otherwise looks at everything
%        (defaults to -7.0)
%  OUTPUTS:
%   -FI: structure with information related to F-I curve
%     -FI.F is list of frequencies in Hz
%     -FI.I is list of currents in nA
%     -FI.Slope is slope of non-zero part of FI curve
%     -FI.Offset is offset of non-zero part of FI curve
%     -FI.Thresh_I is estimated spike threshold current
%     -FI.Analyze is list of spike analysis structures for each trace
%    The remaining are value/voltage pairs for the locations:
%      -PreMaxCurve = maximum curvature before spike
%      -Deriv = maximum derivative
%      -SpikeHeight and SpikeWidth (hopefully obvious)

%defaults:
tCutoff = -7.0;
PlotVar = false;

if(nargin < 1)
  Okay = false;
  BadNumArg = false;
elseif(ischar(varargin{1}))   %specified FileName
  if(nargin < 2 || nargin > 4)
    BadNumArg = true;
  else
    BadNumArg = false;
    FileName = varargin{1};
    CellName = varargin{2};
    if(nargin >= 3 && length(varargin{3}) > 0)
      PlotVar = varargin{3};
    end
    if(nargin >= 4 && length(varargin{4}) > 0)
      tCutoff = varargin{4};
    end
    [t, v, I, Okay] = LoadCleanIntrinsicData(FileName, CellName);
  end
else
  if(nargin < 3 || nargin > 5)
    BadNumArg = true;
  else
    BadNumArg = false;
    Okay = true;
    t = varargin{1};
    v = varargin{2};
    I = varargin{3};
    if(nargin >= 4 && length(varargin{4}) > 0)
      PlotVar = varargin{4};
    end
    if(nargin >= 5 && length(varargin{5}) > 0)
      tCutoff = varargin{5};
    end
  end
end
if(BadNumArg)
  error('Incorrect number of input arguments.  Run "help GetFICurve"')
elseif(~Okay)
  FI.F = [];
  FI.I = [];
  FI.Slope = NaN;
  FI.Offset = NaN;
  FI.Thresh_I = NaN;
  FI.Analyze = [];
  FI.PreMaxCurve.K = [];
  FI.PreMaxCurve.V = [];
  FI.Deriv.D = [];
  FI.Deriv.V = [];
  FI.SpikeHeight.H = [];
  FI.SpikeHeight.V = [];
  FI.SpikeWidth.W = [];
  FI.SpikeWidth.V = [];
  return
end

IndivSpikePlotVar = false;
OverallSpikePlotVar = true;
NoShape = false;
NumTraces = size(v,2);

%First find when the current turns on and off
Ind = find(I(:,end) > .1);
On = Ind(1) + 1;
Off = Ind(end) - 1;
if(tCutoff > 0)
  AvgOn = On;
  AvgOff = On + round(1000 * tCutoff / (t(2) - t(1)));
elseif(tCutoff < 0)
  AvgOff = Off;
  AvgOn = Off + round(1000 * tCutoff / (t(2) - t(1)));
else
  AvgOn = On;
  AvgOff = Off;
end

STimes = cell(NumTraces,1);
SMaxes = cell(NumTraces,1);
FI_F = [];
FI_I = [];
FI_K_K = [];
FI_K_V = [];
FI_D_D = [];
FI_D_V = [];
FI_SH_H = [];
FI_SH_V = [];
FI_SW_W = [];
FI_SW_V = [];
FI_Analyze = [];
for n = 1:NumTraces
  Analyze = AnalyzeWaveform3(t(AvgOn:AvgOff), v(AvgOn:AvgOff,n), ...
			     IndivSpikePlotVar, NoShape);

  STimes{n} = Analyze.Spike.Times;
  SMaxes{n} = Analyze.Spike.MaxV.V.List;
  Num = length(STimes{n});
  if(Num == 0)
    Freq = 0;
  elseif(Num == 1)
    Freq = 0; %Analyze.Spike.Freq;
  else
    Freq = 1000 * (Num - 1) / (STimes{n}(end) - STimes{n}(1));
  end
  FI_F = [FI_F, Freq];
  FI_I = [FI_I, mean(I(On:Off,n))];
  [FI_K_K, FI_K_V] = GetMeanFIProps(FI_K_K, FI_K_V, Analyze, ...
				    {'PreMaxCurve', 'K'});
  [FI_D_D, FI_D_V] = GetMeanFIProps(FI_D_D, FI_D_V, Analyze, ...
				    {'MaxDeriv', 'DV'});
  [FI_SH_H, FI_SH_V] = GetMeanFIProps(FI_SH_H, FI_SH_V, Analyze, ...
				      {'MaxV', 'V'}, ...
				      {'PreMaxCurve', 'V'});
  [FI_SW_W, FI_SW_V] = GetMeanFIProps(FI_SW_W, FI_SW_V, Analyze, ...
				      {'PostMinV', 't'}, ...
				      {'PreMaxCurve', 't'});				       
  FI_Analyze = [FI_Analyze, Analyze];
end

Ind = find(FI_F > 0 & FI_I > 0, 1);
if(length(Ind) == 0)
  Slope = NaN;
  Offset = NaN;
  ThreshCurrent = NaN;
elseif(Ind > length(FI_F) - 2)
  Slope = NaN;
  Offset = NaN;
  ThreshCurrent = NaN;
else
  for n = length(FI_I):-1:(Ind+2)
    FitVals = polyfit(FI_I(Ind:n), FI_F(Ind:n), 1);
    Residuals = FI_F(Ind:n) - (FitVals(1) * FI_I(Ind:n) + FitVals(2));
    ChiSquared = Residuals * Residuals' / (length(Residuals) - 2);
    if(ChiSquared < 2)
      break
    end
  end
  
  Slope = FitVals(1);
  Offset = FitVals(2);
  ThreshCurrent = -Offset / Slope;
end

if(PlotVar)
  Ind = [strfind(FileName, '/'), strfind(FileName, '\')];
  if(length(Ind) > 0)
    FName = FileName((Ind(end)+1):end);
  else
    FName = FileName;
  end
  Ind = strfind(FName, '_');
  for n=length(Ind):-1:1
    m = Ind(n);
    FName = [FName(1:(m-1)), '\', FName(m:end)];
  end

  TitleStr = FName;
  h = NamedFigure(TitleStr);
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(0.001 * t, v)
  if(OverallSpikePlotVar)
    hold on
    for n = 1:NumTraces
      plot(0.001 * STimes{n}, SMaxes{n}, 'k.', 'MarkerSize', 5)
    end
    hold off
  end
  xlabel('Time (s)', 'FontSize', 18)
  ylabel('Voltage (mV)', 'FontSize', 18)
  title(TitleStr, 'FontSize', 18)
  
  TitleStr = ['FI Plot for ', FName];
  h = NamedFigure(TitleStr);
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(FI_I, FI_F, 'rx', 'MarkerSize', 10)
  hold on
  plot(FI_I, Offset + Slope * FI_I, 'b-')
  hold off
  xlabel('Current (nA)', 'FontSize', 18)
  ylabel('Frequency (Hz)', 'FontSize', 18)
  title(TitleStr, 'FontSize', 18)
  
end

%if(Slope < 0)
%  keyboard
%end

FI.F = FI_F;
FI.I = FI_I;
FI.Slope = Slope;
FI.Offset = Offset;
FI.Thresh_I = ThreshCurrent;
FI.Analyze = FI_Analyze;
FI.PreMaxCurve.K = FI_K_K;
FI.PreMaxCurve.V = FI_K_V;
FI.Deriv.D = FI_D_D;
FI.Deriv.V = FI_D_V;
FI.SpikeHeight.H = FI_SH_H;
FI.SpikeHeight.V = FI_SH_V;
FI.SpikeWidth.W = FI_SW_W;
FI.SpikeWidth.V = FI_SW_V;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PropList, VList] = GetMeanFIProps(PropList, VList, ...
					   Analyze, SNames1, SNames2)
if(nargin == 4)
  TempProp = Analyze.Spike.(SNames1{1}).(SNames1{2}).Mean;
  TempV =Analyze.Spike.(SNames1{1}).('V').Mean;
  PropList = [PropList, TempProp];
  VList = [VList, TempV];
else  %it's a difference between two properties
  if(strcmp(SNames1{2}, 't'))   %it's a time-difference
    TempProp1 = Analyze.Spike.(SNames1{1}).(SNames1{2});
    TempProp2 = Analyze.Spike.(SNames2{1}).(SNames2{2});
  else                      %it's not a time difference
    TempProp1 = Analyze.Spike.(SNames1{1}).(SNames1{2}).List;
    TempProp2 = Analyze.Spike.(SNames2{1}).(SNames2{2}).List;
  end
  if(length(TempProp1) == 0)% causes problems with older Matlab versions
    TempProp = NaN;
  else
    TempProp = mean(TempProp1 - TempProp2);
  end
  TempV =Analyze.Spike.PreMaxCurve.V.Mean;
  PropList = [PropList, TempProp];
  VList = [VList, TempV];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, v, I] = GetFIData(FileName)
AbfS = LoadAbf(FileName);

Current = [];
Voltage = [];
FNames = fieldnames(AbfS.Units);
for n = 1:length(FNames)
  Unit = AbfS.Units.(FNames{n});
  if(strcmp(Unit, 'mV'))
    Voltage = [Voltage, n];
  elseif(strcmp(Unit, 'nA'))
    Current = [Current, n];
  end
end
if(length(Current) ~= 1 | length(Voltage) ~= 1)
  AbfS.Units
  error(['Error opening ', FileName, ...
	 ':  Expected one voltage file and one current file.'])
end
v = AbfS.Data.(FNames{Voltage});
I = AbfS.Data.(FNames{Current});
t = AbfS.Time * 1000;  %convert to ms

[t, v, I, DCC_Info] = CleanDCC(t, v, I);
if(~isfinite(DCC_Info.DCC_Freq))
  disp('Warning:  weird DCC signal!')
  DCC_Info
end

SmoothTime = .5;
if(t(2) - t(1) < SmoothTime)
  n = 2;
  while(t(n) - t(1) < SmoothTime)
    n = n + 2;
  end
  [B,A] = butter(2, 2 / n,'low');
  v = filtfilt(B, A, v);
  I = filtfilt(B, A, I);
end
return