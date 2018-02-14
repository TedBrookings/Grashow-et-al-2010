function IV = GetIVCurve(varargin)
% Note that this is poorly-named:  should be GetVICurve.m
% IV = GetIVCurve(FileName, CellName, PlotVar, tCutOff)
%        ---OR---
% IV = GetIVCurve(t, v, I, PlotVar, tCutoff)
% if no parameters are passed, a structure filled with NaNs is returned.
%  INPUTS:
%   -t is time in ms
%   -v is voltage in mV
%   -I is current in nA
%   -FileName is the name of an .abf file with IV data
%   -CellName is the name of the requested cell within the
%      file specified by FileName (e.g. "top" or "bot")
%    OPTIONAL:
%     -PlotVar set to true/1 to plot the result
%     (defaults to false/0)
%     -tCutoff sets time (in s) during current injection to look at spikes
%        if tCutoff > 0, looks at first tCutoff seconds of data
%        if tCutoff < 0, looks at last abs(tCutoff) seconds of data
%        otherwise looks at everything
%        (defaults to 3.0)
%      NOTE:  returned voltage is **MIN**  if tCutoff > 0
%                                 **MEAN** otherwise
%  OUTPUTS:
%   -IV: structure with information related to I-V curve
%     -IV.I is list of currents in nA
%     -IV.V is list of voltages in mV
%     -IV.R is input impedence in megaohms
%     -IV.VRest is resting membrane potential in mV, if the
%       recording has non-spiking data with roughly zero input
%       current.  Otherwise it is NaN
%     -IV.VIntercept is estimated resting membrane potential in
%       mV, as linearly fit from I-V curve
%     -IV.ChiSquared

%defaults:
tCutoff = 3.0;
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
  error('Incorrect number of input arguments.  Run "help GetIVCurve"')
elseif(~Okay)
  IV.I = NaN;
  IV.V = NaN;
  IV.R = NaN;
  IV.VRest = NaN;
  IV.VIntercept = NaN;
  IV.ChiSquared = NaN; 
  return
end

%First find when the current turns on and off
BaseI = mean(I(2,:));

ABS_I = abs(I - BaseI);
[Max1, Ind1] = max(ABS_I);
[MaxI, Ind2] = max(Max1);

Ind = find(ABS_I(:,Ind2) > .15);
if(length(Ind) < 100)
  NumTraces = size(I,2);
  %I has a big blip in it that is screwing things up in this trace
  ABS_I = abs(I(:,[1:(Ind2-1),(Ind2+1):NumTraces]));
  [Max1, Ind1] = max(ABS_I);
  [MaxI, Ind2] = max(Max1);

  Ind = find(ABS_I(:,Ind2) > .15);
end
On = Ind(1) + 1;
Off = Ind(end) - 1;
if(tCutoff > 0)
  AvgOn = On;
  AvgOff = On + round(1000 * tCutoff / (t(2) - t(1)));
  vFunc = @mean;
elseif(tCutoff < 0)
  AvgOff = Off;
  AvgOn = Off + round(1000 * tCutoff / (t(2) - t(1)));
  vFunc = @min;
else
  AvgOn = On;
  AvgOff = Off;
  vFunc = @mean;
end

ShowPlot = false;
NoShape = true;
%Get the resting membrane potential
IMean = abs(mean(I(On:Off, :)));
[IMin, RestInd] = min(IMean);
if(IMin > .05)
  V0 = NaN;  %Nothing is close to zero current
  SpikesAtRest = true;  %Just to be sure, check for spikes
else
  %Check if the cell is spiking at rest
  Analyze = AnalyzeWaveform3(t, v(:,RestInd), ShowPlot, NoShape);
  SpikesAtRest = (length(Analyze.Spike.Times) > 0);
  if(SpikesAtRest)
    V0 = NaN;
  else
    V0 = mean(v(:,RestInd));
  end
end

IV_I = [];
IV_V = [];

%Loop through the traces
if(SpikesAtRest)
  Stop_n = RestInd - 1;
else
  Stop_n = RestInd;
end
for n = 1:Stop_n
  TempI = mean(I(On:Off,n));
  if(TempI > .05)
    break
  end
  if(SpikesAtRest)
    Mid = round(.5*(On + Off));

    Analyze = AnalyzeWaveform3(t(Mid:Off), v(Mid:Off,n), ShowPlot, NoShape);
    if(length(Analyze.Spike.Times) > 0)
      %Analyze = AnalyzeWaveform2(t(Mid:Off), v(Mid:Off,n), true, NoShape);
      %keyboard
      break
    end
  end
  
  IV_I = [IV_I, TempI];
  IV_V = [IV_V, vFunc(v(AvgOn:AvgOff,n))];
end

if(length(IV_I) < 3)
  disp('I-V curve has fewer than three points.  Returning "empty" struct.')
  IV.I = NaN;
  IV.V = NaN;
  IV.R = NaN;
  IV.VRest = NaN;
  IV.VIntercept = NaN;
  IV.ChiSquared = NaN; 
  return
end
SaveLineVals = polyfit(IV_I, IV_V, 1);

Residuals = IV_V - (SaveLineVals(1) * IV_I + SaveLineVals(2));
SaveChiSquared = Residuals * Residuals' / (length(Residuals) - 2);
if(SaveChiSquared > 1)
  for n = (length(IV_I)-1):3
    LineVals = polyfit(IV_I(1:n), IV_V(1:n), 1);
    Residuals = IV_V(1:n) - (LineVals(1) * IV_I(1:n) + LineVals(2));
    ChiSquared = Residuals * Residuals' / (length(Residuals) - 2);
    
    if(ChiSquared < SaveChiSquared)
      SaveChiSquared = ChiSquared;
      SaveLineVals = LineVals;
    end
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
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(0.001 * t, v)
  hold on
  plot(0.001 * [t(1), t(end)], [IV_V; IV_V]);
  hold off
  xlabel('Time (s)', 'FontSize', 18)
  ylabel('Voltage (mV)', 'FontSize', 18)
  title(TitleStr, 'FontSize', 18)
  
  Title2 = [TitleStr, ': Injected Current'];
  h = NamedFigure(Title2);
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(0.001 * t, I)
  xlabel('Time (s)', 'FontSize', 18)
  ylabel('Current (nA)', 'FontSize', 18)
  title(Title2, 'FontSize', 18)
  
  Title3 = [TitleStr, ': IV Curve and Fit'];
  h = NamedFigure(Title3);
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(IV_I, IV_V, 'bx', 'MarkerSize', 10)
  hold on
  plot(IV_I, SaveLineVals(1) * IV_I + SaveLineVals(2), 'r-');
  hold off
  xlabel('Current (nA)', 'FontSize', 18)
  ylabel('Voltage (mV)', 'FontSize', 18)
  title(Title3, 'FontSize', 18)
end

IV.I = IV_I;
IV.V = IV_V;
IV.R = SaveLineVals(1);
IV.VRest = V0;
IV.VIntercept = SaveLineVals(2);
IV.ChiSquared = SaveChiSquared;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = GetIVAndTau(t, I, v, V0, On, Off, Trace)
IV_I = mean(I(On:Off,Trace));

xData = t(On:Off);
yData = v(On:Off, Trace);

DeltaV = yData(end) - V0;
if(DeltaV < 0)
  VRange1 = [1.5 * DeltaV, .7 * DeltaV];
  VRange2 = [.9 * DeltaV, -.5 * DeltaV];
else
  VRange1 = [.7 * DeltaV, 1.5 * DeltaV];
  VRange2 = [-.5 * DeltaV, .9 * DeltaV];
end

fHandle = @VCurve;
StartParams = [V0, .8 * DeltaV, .005, .2 * DeltaV, .1];
ParamRanges = [[V0, V0]; VRange1; [.001, .3]; VRange2; [.02, 2]];

DeltaT = t(2) - t(1);
Tau1 = DeltaT / log(abs((yData(3) - yData(1))/(yData(2) - yData(1))));
DeltaV1 = yData(2) - yData(1);
DeltaV2 = DeltaV - DeltaV1;
Tau2 = 10 * Tau1;
%StartParams = [V0, DeltaV1, Tau1, DeltaV2, Tau2];

VRange1 = sort([-DeltaV * .5, DeltaV * 1.5]);
VRange2 = sort([-DeltaV * .5, DeltaV * 1.5]);
TauRange1 = [Tau1 / 10, Tau1 * 10];
TauRange2 = [Tau2 / 10, Tau2 * 10];
%ParamRanges = [[V0, V0]; VRange1; TauRange1; VRange2; TauRange2];


Tol = 1e-4;
fTol = 1e-5;
DerivTol = 1e-7;

Verbose = false;

try
  Params = FitChiSquared(fHandle, StartParams, ParamRanges, ...
			 xData, yData, ...
			 Tol, fTol, DerivTol, ...
			 Verbose);
catch
  
  disp('Error fitting in GetIVCurve.m, GetIVAndTau()')
  keyboard
end
DebugFit = true;
if(DebugFit)
  Params
  hold off
  plot(xData, yData)
  hold on
  yFit = VCurve(Params, xData);
  plot(xData, yFit, 'r-');
  
  keyboard
end

  
DeltaV1 = Params(2);
Tau1 = Params(3);
DeltaV2 = Params(4);
Tau2 = Params(5);

IV_V = V0 + DeltaV1 + DeltaV2;

varargout = {IV_I, IV_V, DeltaV1, Tau1, DeltaV2, Tau2};
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V, varargout] = VCurve(Parameters, t)
%Parameters(1) = V0  (fixed by data)
%Parameters(2) = DeltaV1
%Parameters(3) = Tau1
%Parameters(4) = DeltaV2
%Parameters(5) = Tau2

Temp1 = 1 - exp(-(t - t(1)) / Parameters(3));
Temp2 = 1 - exp(-(t - t(1)) / Parameters(5));
V = Parameters(1) + Parameters(2) * Temp1 ...
    + Parameters(4) * Temp2;
return
if(nargout > 1)
  Deriv = [ones(size(t)), Temp1, ...
	   Parameters(2) / Parameters(3)^2 * (Temp1 .* (t - t(1))), ...
	   Temp2, ...
	   Parameters(4) / Parameters(5)^2 * (Temp2 .* (t - t(1)))];
  varargout = {Deriv};
end
return