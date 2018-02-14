function Ramp = GetRampData(varargin)
% Ramp = GetRampData(FileName, CellName, PlotVar, tCutoff)
%        ---OR---
% Ramp = GetRampData(t, v, I, PlotVar, tCutoff)
% if no parameters are passed, a structure filled with NaNs is returned.
% Usually, you want Ramp.PreMaxCurve.V (threshold voltage)
%                or Ramp.PreMaxCurve.I (threshold current)
%  INPUT PARAMETERS:
%   -t is time in ms
%   -v is voltage in mV
%   -I is current in nA
%   -FileName is the name of an .abf file with IV data
%   -CellName is the name of the requested cell within the
%      file specified by FileName (e.g. "top" or "bot")
%    OPTIONAL:
%     -PlotVar set to 1/true to plot the result
%     (defaults to 0/false)
%     -tCutoff sets time (in s) during current injection to look at spikes
%        if tCutoff > 0, looks at first tCutoff seconds of data
%        if tCutoff < 0, looks at last abs(tCutoff) seconds of data
%        otherwise looks at everything
%        (defaults to 0.030)
%  OUTPUT PARAMETERS:
%   -Ramp: structure with information related to current ramp data
%     -Ramp.Analyze:  structure produced by AnalyzeWaveform3.m
%     -Ramp.SpikeHeight:  voltage of spike peak, minus voltage at
%          point of maximum curvature before spike
%     -Ramp.SpikeWidth:  time of spike peak, minus time at point of
%          maximum curvature before spike
%     -Ramp.TenMs: structure with I/V data for point 10 ms before spike peak
%       -Ramp.TenMs.V:  voltage in mV
%       -Ramp.TenMs.I:  current in nA
%     There are several other structures, with hopefully
%      self-descriptive names.  They all contain, a .t, .V and .I
%      field, corresponding to the time that they are sampled, the
%      voltage, and the current.  They may also have a fourth field
%      (e.g. K for maximum curvature).  They are:
%        -Ramp.MaxV, Ramp.MaxDeriv, Ramp.MinDeriv,
%         Ramp.PreMaxCurve, Ramp.PostMaxCurve

Ramp.MaxV.V = NaN;
Ramp.MaxV.t = NaN;
Ramp.MaxV.I = NaN;
Ramp.MaxDeriv.V = NaN;
Ramp.MaxDeriv.DV = NaN;
Ramp.MaxDeriv.t = NaN;
Ramp.MinDeriv.V = NaN;
Ramp.MinDeriv.I = NaN;
Ramp.MinDeriv.DV = NaN;
Ramp.MinDeriv.t = NaN;
Ramp.MinDeriv.I = NaN;
Ramp.PreMaxCurve.t = NaN;
Ramp.PreMaxCurve.V = NaN;
Ramp.PreMaxCurve.K = NaN;
Ramp.PreMaxCurve.I = NaN;
Ramp.PostMaxCurve.t = NaN;
Ramp.PostMaxCurve.V = NaN;
Ramp.PostMaxCurve.K = NaN;
Ramp.PostMaxCurve.I = NaN;

%defaults:
tCutoff = 0.030;
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
  error('Incorrect number of input arguments.  Run "help GetRampData"')
elseif(~Okay)
  Ramp.SpikeHeight = NaN;
  Ramp.SpikeWidth = NaN;
  Ramp.TenMs.I = NaN;
  Ramp.TenMs.V = NaN;
  Ramp.Analyze = NaN;
  return
end

NoShape = false;
FirstOnly = true;

%First find when the current turns on and off
BaseI = mean(I(2:10));
Ind = find(abs(I - BaseI) > 0.1, 1);
if(tCutoff > 0)
  On = Ind(1) + round(1000 * tCutoff / (t(2) - t(1)));
elseif(tCutoff < 0)
  On = length(t) + round(1000 * tCutoff / (t(2) - t(1)));
else
  On = Ind(1);
end

numSweeps = size(v, 2);
if(numSweeps == 1)
  Ramp = AnalyzeSweep(t(On:end), v(On:end), I(On:end), Ramp);
else
  Ramp = AnalyzeSweep(t(On:end), v(On:end,1), I(On:end), Ramp);
  for n = 2:numSweeps
    Ramp = [Ramp, AnalyzeSweep(t(On:end), v(On:end,n), I(On:end), ...
			       Ramp(end))];
  end
  
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
  whitebg(h, 'k')
  %set(h, 'WindowStyle', 'docked');
  hold off
  plot(0.001 * t, v)
  hold on
  for n = 1:numSweeps
    plotSpike(Ramp(n))
  end
  hold off
  xlabel('Time (s)', 'FontSize', 18)
  ylabel('Voltage (mV)', 'FontSize', 18)
  title(TitleStr, 'FontSize', 18)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ramp = AnalyzeSweep(t, v, I, Ramp, PlotVar)
if(nargin < 5)
  PlotVar = false;
end
NoShape = false;
FirstOnly = true;

Analyze = AnalyzeWaveform3(t, v, PlotVar, NoShape, FirstOnly);
FNames = fieldnames(Ramp);
for n = 1:length(FNames)
  FName = FNames{n};
  try
    FirstSpikeT = Analyze.Spike.(FName).t;
  catch
    continue
  end
  if(length(FirstSpikeT) ~= 1)
    AnalyzeWaveform3(t, v, true, NoShape, FirstOnly);
    FirstSpikeT
    error('Incorrect number of spikes found')
  end
  [MinDiff, ind] = min(abs(t - FirstSpikeT));
  Ramp.(FName).I = I(ind);
  
  SubFNames = fieldnames(Ramp.(FName));
  for m = 1:length(SubFNames)
    SubFName = SubFNames{m};
    if(strcmp(SubFName, 'I'))
      continue
    elseif(strcmp(SubFName, 't'))
        Ramp.(FName).(SubFName) = Analyze.Spike.(FName).(SubFName);
    else
        Ramp.(FName).(SubFName) = Analyze.Spike.(FName).(SubFName).List;
    end
  end
end
Ramp.SpikeHeight = Ramp.MaxV.V - Ramp.PreMaxCurve.V;
Ramp.SpikeWidth = Analyze.Spike.PostMinV.t - Ramp.PreMaxCurve.t;

tTenMs = Ramp.MaxV.t - 10;
[MinVal, Ind] = min(abs(t - tTenMs));
%[MinVal, Ind] = min(abs(t - Ramp.PreMaxCurve.t));
Ramp.TenMs.I = I(Ind);
Ramp.TenMs.V = v(Ind);

Ramp.Analyze = Analyze;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSpike(Ramp)
spikeT = 0.001 * Ramp.MaxV.t;
maxV = Ramp.MaxV.V;
yRange = maxV + Ramp.SpikeHeight * [-2, 1];
plot([spikeT, spikeT], yRange, 'r-')
plot(0.001 * Ramp.PreMaxCurve.t, Ramp.PreMaxCurve.V, 'wx', ...
     'MarkerSize', 7)

return