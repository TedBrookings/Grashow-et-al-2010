function StructOut = AnalyzeWaveform(t, v, varargin)
% StructOut = AnalyzeWaveform(t, v, PlotSubject, NoShape, FirstOnly)
% Analyzes a single voltage waveform, looking for spikes
%    and bursts, and calculating relevant frequencies.
%
%  INPUT PARAMETERS:
%   -t is time in ms
%   -V is voltage in mV
%    OPTIONAL:
%     -PlotSubject should be set to true[false] to produce[suppress]
%       plots of waveforms/analysis.  Alternatively, it can be set
%       to a string to aid it titling plots (e.g. 'Exp #71')
%       PlotVar defaults to false
%     -NoShape: set to true[false] to suppress[enforce] calculation
%      of spike shape.  Defaults to false.
%     -FirstOnly: set to true[false] to stop[continue] analysis
%      after the first spike is detected.  Defaults to false.
%  OUTPUT PARAMETERS:
%   -StructOut.Spike:  structure with spike information
%      -Spike.Freq is overall spiking frequency (in Hz)
%      -Spike.Times is a plain list of spike times (in ms)
%      -Spike.Intervals is a list of interspike intervals (in ms)
%      -Spike.Frequencies is a list of instantaneous frequencies (in Hz)
%      Shape information structures (should be self-descriptive)
%      -Spike.MaxV, Spike.MaxDeriv, Spike.MinDeriv, Spike.PreMinV,
%       Spike.PostMinV, Spike.PreMaxCurve, Spike.PostMaxCurve
%           Each contains a list of times/voltage points, and if relevant
%           another quantity (such as K for curvatures)
%   -StructOut.SlowWave:  structure with slow-wave information
%      -SlowWave.Freq: frequency of the dominant slow-wave
%       component (in Hz)
%      -SlowWave.Sigma: (very crude) measure of the importance
%       of the slow-wave frequency in the power spectrum
%      -SlowWave.Corr: autocorrelation at slow-wave period
%      -SlowWave.Spectrum: structure with spectrum information
%        -Spectrum.Freq:  list of analyzed frequencies
%        -Spectrum.Power: length NumFreq list of average powers of
%           waveform with spikes removed.
%   -StructOut.Burst:  structure with burst information
%      -Burst.Freq is burst frequency (in Hz)
%      -Burst.SpikeFreq is within-burst spike frequency (in Hz)
%      -Burst.DutyCycle is the average burst duration/period
%      -Burst.Times is a plain list of burst times (in ms)
%      -Burst.Durations is a list of burst durations (in ms)
%      -Burst.NumSpikes is a list of spikes per burst
%      -Burst.SpikeFrequencies is a list of spike frequencies (in Hz)
%      -Burst.InterBurstIntervals is a list of inter-burst
%       intervals (in ms)
%   -StructOut.MedianV:  the median of the voltage trace.  If the
%      cell is silent, it should be the resting potential,
%      otherwise, who knows...
%
%List structures usually will have a Name.List element, as well as
%  Name.Mean, Name.StdDev, Name.Variance, Name.CoefOfVar
%  (a few are just plain lists)
%If a feature is not detected, relevant frequencies are set to
%  zero, and relevant lists are empty
%
%NOTE for future:  would benefit enormously by changing to .mex
tic

if(nargin < 3)
  PlotSubject = false;
  NoShape = false;
  FirstOnly = false;
else
  PlotSubject = varargin{1};
  if(nargin < 4)
    NoShape = false;
    FirstOnly = false;
  else
    NoShape = varargin{2};
    if(nargin < 5)
      FirstOnly = false;
    else
      FirstOnly = varargin{3};
    end
  end
end

%First get the spike times
Spike = GetSpikeTimes(t, v, PlotSubject, NoShape, FirstOnly);

%Next get the overall spike frequency
Spike.Freq = GetSpikeFrequency(Spike.Times, t);

if(FirstOnly)
  SlowWave = [];
  Burst = [];
else
  %keyboard
  SlowWave = AnalyzeSlowWave(t, v, Spike, PlotSubject);
  Burst = AnalyzeBurst(Spike, SlowWave, t);
  SlowWave.Phases = [];  %Reduce storage demand
end

%Structify (add info about mean, variance, etc) various lists
Spike.Intervals = StructifyList(Spike.Intervals);
Spike.Frequencies = StructifyList(Spike.Frequencies);
%Maybe should structify other shape information
if(~NoShape)
  Spike.MaxV.V = StructifyList(Spike.MaxV.V);
  Spike.MaxDeriv.V = StructifyList(Spike.MaxDeriv.V);
  Spike.MaxDeriv.DV = StructifyList(Spike.MaxDeriv.DV);
  Spike.MinDeriv.V = StructifyList(Spike.MinDeriv.V);
  Spike.MinDeriv.DV = StructifyList(Spike.MinDeriv.DV);
  Spike.PreMinV.V = StructifyList(Spike.PreMinV.V);
  Spike.PostMinV.V = StructifyList(Spike.PostMinV.V);
  Spike.PreMaxCurve.V = StructifyList(Spike.PreMaxCurve.V);
  Spike.PreMaxCurve.K = StructifyList(Spike.PreMaxCurve.K);
  Spike.PostMaxCurve.V = StructifyList(Spike.PostMaxCurve.V);
  Spike.PostMaxCurve.K = StructifyList(Spike.PostMaxCurve.K);
end

Burst.Durations = StructifyList(Burst.Durations);
Burst.SpikesPerBurst = StructifyList(Burst.SpikesPerBurst);
Burst.InterBurstIntervals = StructifyList(Burst.InterBurstIntervals);
Burst.SpikeFrequencies = StructifyList(Burst.SpikeFrequencies);

StructOut.Spike = Spike;
StructOut.SlowWave = SlowWave;
StructOut.Burst = Burst;
StructOut.MedianV = median(v);

if(DoPlot(PlotSubject))
  PlotGetSpikes(t, v, Spike, Burst, PlotSubject);
end
Toc
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Spike = GetSpikeTimes(t, V, varargin)
%Find spikes as traces that travel clockwise around a fixed point
% on a graph of V vs dV/dt

DebugPlots = false;
%Determine whether to plot and/or calculate spike-shape data
if(nargin < 3)
  PlotSubject = false;
  NoShape = false;
  FirstOnly = false;
else
  PlotSubject = varargin{1};
  if(nargin < 4)
    NoShape = false;
    FirstOnly = false;
  else
    NoShape = varargin{2};
    if(nargin < 5)
      FirstOnly = false;
    else
      FirstOnly = varargin{3};
    end
  end
end

%Calculate derivatives as needed
DeltaT = t(2) - t(1);
if(NoShape)
  Deriv = PolyDeriv(V, DeltaT);
  Deriv2 = 0;
else
  [Deriv, Deriv2] = PolyDeriv(V, DeltaT);
end

%High-pass filter voltage for this analysis
OldV = V;
SmoothTime = 30;  %ms
if(t(2) - t(1) < SmoothTime)
  n = 2;
  while(t(n) - t(1) < SmoothTime)
    n = n + 2;
  end
  [B,A] = butter(2, 2 / n,'high');
  V = filtfilt(B, A, OldV);
end

[VMid, dVMid] = GetInteriorPoint(V, Deriv);

[n1List, n2List, BadCount] = GetIndLists(V, Deriv, VMid, dVMid, ...
					 NoShape, FirstOnly);

[n1List2, n2List2, BadCount2, VMid2, dVMid2] ...
    = RecurseGetSpikes(V, Deriv, n1List, n2List, BadCount, NoShape, FirstOnly);
if(BadCount > .2 * length(n1List))
  if(BadCount2 > .2 * length(n1List2))
    n1List = [];
    n2List = [];
  else
    n1List = n1List2;
    n2List = n2List2;
    VMid = VMid2;
    dVMid = dVMid2;
  end
elseif(BadCount2 <= .2 * length(n1List2))
  if(BadCount2 / length(n1List2) < BadCount / length(n1List))
    n1List = n1List2;
    n2List = n2List2;
    VMid = VMid2;
    dVMid = dVMid2;
  end
end

%now spikes are all bracketed between n1List and n2List

%Get spike shape (either dummy struct, or filled with data)
Spike = GetSpikeShape(n1List, n2List, t, OldV, Deriv, Deriv2, NoShape);

if(DoPlot(PlotSubject) & DebugPlots)
  PlotGetSpikeTimes(t, V, Deriv, PlotSubject);
end

if(length(Spike.Times) > 0)
  Spike.Intervals = Spike.Times(2:end) - Spike.Times(1:(end-1));
  Spike.Frequencies = 1000 ./ Spike.Intervals;
else
  Spike.Intervals = [];
  Spike.Frequencies = [];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second try, maybe a few weird points screwed-up spike detection
function [n1List, n2List, BadCount, VMid, dVMid] ...
    = RecurseGetSpikes(V, Deriv, n1List, n2List, BadCount, NoShape, FirstOnly)
NumV = length(V);
Cutoff = .1 * NumV;
if(length(n1List) > Cutoff | BadCount > Cutoff)
  n1List = [];
  n2List = [];
  VMid = Inf;
  dVMid = Inf;
  BadCount = Inf;
  return
end
GoodInd = 1:NumV;
for Ind = length(n1List):-1:1
  n1 = n1List(Ind) - 1;
  n2 = n2List(Ind) + 1;
  GoodInd = [GoodInd(1:n1), GoodInd(n2:end)];
end

[VMid, dVMid] = GetInteriorPoint(V(GoodInd), Deriv(GoodInd));

[n1List, n2List, BadCount] = GetIndLists(V, Deriv, VMid, dVMid, ...
					 NoShape, FirstOnly);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n1List, n2List, BadCount] = GetIndLists(V, Deriv, VMid, dVMid, ...
						  NoShape, FirstOnly)
%Loop through data, looking for spikes
NumV = length(V);
n1List = [];
n2List = [];
n = 3;
BadCount = 0;  %BadCount keeps track of the number of times V gets
               %large without ever achieving a large derivative.
               %If BadCount is large, it means there are no spikes,
               %just random voltage fluctuations
while(n <= NumV)
  if(V(n) < VMid)
    n = n + 1;
  else
    if(Deriv(n) > dVMid)
      n1 = n - 1;
    elseif(dVMid < Deriv(n-1) + (VMid - V(n-1)) * ...
				 (Deriv(n) - Deriv(n-1)) / (V(n) - V(n-1)))
      n1 = n - 2;
    else  %Probably on a bad trajectory
      n1 = n - 1;
      n = n + 1;
      if(n >= NumV)
        BadCount = BadCount + 1;
        break;
      end
      while(V(n) >= VMid & Deriv(n) <= dVMid)
	n = n + 1;
	if(n >= NumV)
	  BadCount = BadCount + 1;
	  break
	end
      end
      if(n >= NumV)
	break
      end
      if(V(n) < VMid)
	BadCount = BadCount + 1;
	n = n + 1;
	continue
      end
    end

    n2 = n + 1;

    while(n2 <= NumV)
      if(V(n2) > VMid)
	n2 = n2 + 1;
      else
	break
      end
    end
    if(n2 > NumV)
      break
    end
    
    %We've bracketed a spike between n1 and n2
    if(NoShape)
      n1List = [n1List, n1];
      n2List = [n2List, n2];
      if(FirstOnly)
	break
      end
      n = n2;
      continue
    end
    
    %We want to get some spike shape info, so extend out n1 and n2
    % to encompass time before and after spike
    while(Deriv(n1) > 0 | Deriv(n1 + 1) > 0 | Deriv(n1 + 2) > 0)
      if(n1 == 1)
	break;
      else
	n1 = n1 - 1;
      end
    end
    n1List = [n1List, n1];
    
    while(Deriv(n2) < 0 | Deriv(n2 - 1) < 0 | Deriv(n2 - 2) < 0)
      if(n2 == NumV)
	break;
      else
	n2 = n2 + 1;
      end
    end
    n = n2;
    n2List = [n2List, n2];
    if(FirstOnly)
      break
    end
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old method, included as a backup.  Finds spikes by looking for
% points where derivative is large (positive) followed quickly by a
% large (negative) derivative.
function Spike = GetSpikeTimesBackup(t, V, varargin)
if(nargin < 3)
  PlotSubject = false;
  NoShape = false;
  FirstOnly = false;
else
  PlotSubject = varargin{1};
  if(nargin < 4)
    NoShape = false;
    FirstOnly = false;
  else
    NoShape = varargin{2};
    if(nargin < 5)
      FirstOnly = false;
    else
      FirstOnly = varargin{3};
    end
  end
end
if(nargin < 2)
  NumV = 0;
else
  NumV = length(V);
end
if(NumV == 0)
  Spike = GetSpikeShape([], [], [], [], [], NoShape);
  Spike.Intervals = [];
  Spike.Frequencies = [];
  return
end

MaxTimeWidth = 30;  %ms
DeltaT = t(2) - t(1);
if(DeltaT < .005)
  MaxTimeWidth = MaxTimeWidth / 1000;
end

if(NoShape)
  Deriv = PolyDeriv(V, DeltaT);
  Deriv2 = 0;
else
  [Deriv, Deriv2] = PolyDeriv(V, DeltaT);
end

AutoCutoffs = true;
if(AutoCutoffs)
  [LowCutoff, HighCutoff] = GetAutoCutoffs(V, Deriv);
else
  LowCutoff = -2;
  HighCutoff = 2;
end

n1List = [];
n2List = [];
n = 2;

if(HighCutoff - LowCutoff < .1)  %no spikes so skip detection
  n = NumV + 1;
end

while(n <= NumV)
  if(Deriv(n) < HighCutoff)
    n = n + 1;
  else  %Found potential beginning of a spike, try to bracket a spike
    n1 = n - 1;
    n = n + 1;
    if(Deriv(n) < HighCutoff)
      n = n + 1;
      continue;
    end
    n2 = n + 1;
    BracketSuccess = false;
    while(n2 <= NumV & t(n2) - t(n1) < MaxTimeWidth)
      if(Deriv(n2) > LowCutoff)
	if(Deriv(n2) > HighCutoff)  %If we've gone back above cutoff,
	  n1 = n2 - 1;                   %reset n1
	end
	n2 = n2 + 1;
      else
	BracketSuccess = true;
	break
      end
    end
    if(~BracketSuccess)
      n = n2 + 1;
      continue;
    end

    %We've bracketed a spike between n1 and n2
    if(NoShape)
      n1List = [n1List, n1];
      n2List = [n2List, n2];
      if(FirstOnly)
	break
      end
      n = n2 + 1;
      continue
    end
    
    %We want to get some spike shape info, so extend n1 and n2
    %until we cross Deriv = 0
    while(n1 >= 1)
      if(Deriv(n1) <= 0)
	break
      else
	n1 = n1 - 1;
      end
    end
    n1 = n1 + 1;
    n1List = [n1List, n1];
    
    while(n2 <= NumV)
      if(Deriv(n2) >= 0)
	break
      else
	n2 = n2 + 1;
      end
    end
    n = n2;
    n2 = n2 - 1;
    n2List = [n2List, n2];
    if(FirstOnly)
      break
    end
  end
end

%now spikes are all bracketed between n1List and n2List

Spike = GetSpikeShape(n1List, n2List, t, V, Deriv, Deriv2, NoShape);

if(DoPlot(PlotSubject))
  PlotGetSpikeTimes(t, V, Deriv, PlotSubject, LowCutoff, HighCutoff);
end

if(length(Spike.Times) > 0)
  Spike.Intervals = Spike.Times(2:end) - Spike.Times(1:(end-1));
  Spike.Frequencies = 1000 ./ Spike.Intervals;
else
  Spike.Intervals = [];
  Spike.Frequencies = [];
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LowCutoff, HighCutoff] = GetAutoCutoffs(V, Deriv)
SortD = sort(Deriv);

Num = round(.01 * length(Deriv));

LowCutoff = 2.0 * SortD(1 + Num);
HighCutoff = 2.0 * SortD(end - Num);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Spike = GetSpikeShape(n1List, n2List, t, V, Deriv, Deriv2, ...
			       NoShape)
NumSpikes = length(n1List);
Spike.n1List = n1List;
Spike.n2List = n2List;

Spike.Times = zeros(1, NumSpikes); %maybe get rid of this

if(~NoShape)
  Spike.MaxV.V = zeros(1, NumSpikes);
  Spike.MaxV.t = zeros(1, NumSpikes);
  Spike.MaxDeriv.V = zeros(1, NumSpikes);
  Spike.MaxDeriv.DV = zeros(1, NumSpikes);
  Spike.MaxDeriv.t = zeros(1, NumSpikes);
  Spike.MinDeriv.V = zeros(1, NumSpikes);
  Spike.MinDeriv.DV = zeros(1, NumSpikes);
  Spike.MinDeriv.t = zeros(1, NumSpikes);
  Spike.PreMinV.V = zeros(1, NumSpikes);
  Spike.PreMinV.t = zeros(1, NumSpikes);
  Spike.PostMinV.V = zeros(1, NumSpikes);
  Spike.PostMinV.t = zeros(1, NumSpikes);
  Spike.PreMaxCurve.t = zeros(1, NumSpikes);
  Spike.PreMaxCurve.V = zeros(1, NumSpikes);
  Spike.PreMaxCurve.K = zeros(1, NumSpikes);
  Spike.PostMaxCurve.t = zeros(1, NumSpikes);
  Spike.PostMaxCurve.V = zeros(1, NumSpikes);
  Spike.PostMaxCurve.K = zeros(1, NumSpikes);
  if(NumSpikes == 0)
    return
  end
  
  K = Deriv2 .* (1 + Deriv.^2).^-1.5;
  K(1:2) = 0;
  K((end-1):end) = 0;
else
  if(NumSpikes == 0)
    return
  end
end
  
%I think get rid of traces
%Shape.t = {};  %trace t
%Shape.V = {};  %trace V
%Shape.DV = {}; %trace dV/dt

for m = 1:NumSpikes
  n1 = n1List(m);
  n2 = n2List(m);

  %Find the moment and voltage of maximum depolarization
  [MaxV, tMaxV, nMaxV] = GetExtremum(V, t, n1, n2);
  if(NoShape)
    Spike.Times(m) = tMaxV;
    continue;
  end
  %Find the max derivative
  [MaxDV, tMaxDV, nMaxDV] = GetExtremum(Deriv, t, n1, n2);
  VMaxDV = interp1(t, V, tMaxDV);
  %Find the min derivative
  [MinDV, tMinDV, nMinDV] = GetExtremum(Deriv, t, n1, n2, 'min', true);
  VMinDV = interp1(t, V, tMinDV);
  
  %Find the max curvature
  [PreMaxK, tPreMaxK] = GetExtremum(K, t, n1, nMaxV, 'max', true);
  VPreMaxK = interp1(t, V, tPreMaxK);
  [PostMaxK, tPostMaxK] = GetExtremum(K, t, nMaxV, n2, 'max', true);
  VPostMaxK = interp1(t, V, tPostMaxK);
  
  %Find minimum voltage before and after spike
  if(n1 > 1)
    n1 = n1 - 1;
  end
  if(n2 < length(V))
    n2 = n2 + 1;
  end
  [PreMinV, tPreMin] = GetExtremum(V, t, n1, n1+3, 'min', true);
  [PostMinV, tPostMin] = GetExtremum(V, t, n2-3, n2, 'min', true);
  
  Spike.Times(m) = tMaxV;

  Spike.MaxV.V(m) = MaxV;
  Spike.MaxV.t(m) = tMaxV;
  Spike.MaxDeriv.V(m) = VMaxDV;
  Spike.MaxDeriv.DV(m) = MaxDV;
  Spike.MaxDeriv.t(m) = tMaxDV;
  Spike.MinDeriv.V(m) = VMinDV;
  Spike.MinDeriv.DV(m) = MinDV;
  Spike.MinDeriv.t(m) = tMinDV;
  Spike.PreMinV.V(m) = PreMinV;
  Spike.PreMinV.t(m) = tPreMin;
  Spike.PostMinV.V(m) = PostMinV;
  Spike.PostMinV.t(m) = tPostMin;
  Spike.PreMaxCurve.V(m) = VPreMaxK;
  Spike.PreMaxCurve.K(m) = PreMaxK;
  Spike.PreMaxCurve.t(m) = tPreMaxK;
  Spike.PostMaxCurve.V(m) = VPostMaxK;
  Spike.PostMaxCurve.K(m) = PostMaxK;
  Spike.PostMaxCurve.t(m) = tPostMaxK;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TWanted = GetTimeAtVoltage(V, T, VWanted)
if(length(V) == 1)
  TWanted = T(1);
  return
end
[MinVal, nBest] = min(abs(V - VWanted));
if(nBest == 1)
  TWanted = T(nBest) + (VWanted - V(nBest)) / (V(nBest + 1) - V(nBest)) ...
	    * (T(nBest + 1) - T(nBest));
elseif(nBest == length(V))
  TWanted = T(nBest) + (VWanted - V(nBest)) / (V(nBest - 1) - V(nBest)) ...
	    * (T(nBest - 1) - T(nBest));  
elseif(sign(V(nBest + 1) - V(nBest)) == sign(VWanted - V(nBest)))
  TWanted = T(nBest) + (VWanted - V(nBest)) / (V(nBest + 1) - V(nBest)) ...
	    * (T(nBest + 1) - T(nBest));
else
  TWanted = T(nBest) + (VWanted - V(nBest)) / (V(nBest - 1) - V(nBest)) ...
	    * (T(nBest - 1) - T(nBest));
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxV, tMax, nMax] = GetExtremum(V, t, n1, n2, MinStr, Simple)
if(nargin < 6)
  Simple = false;
end
if(nargin < 5)
  MinStr = 'max';
end
if(strcmp(lower(MinStr), 'min'))
  [MaxV, nMax] = min(V(n1:n2));
else
  [MaxV, nMax] = max(V(n1:n2));
end
nMax = nMax + n1 - 1;
if(nMax == 1 || nMax == length(t) || Simple)
  tMax = t(nMax);
  return
end

%Refine by modeling trace as parabola
n1 = nMax - 1;
n2 = nMax;
n3 = nMax + 1;

if(V(n1) == V(n2))
  if(V(n2) == V(n3))
    MaxV = V(n2);
    tMax = t(n2);
    return
  else
    tMax = (t(n1) + t(n2))/2;
    Coeff = (V(n2) - V(n3)) / ((t(n2) - tMax)^2 - (t(n3) - tMax)^2);
  end
elseif(V(n2) == V(n3))
  tMax = (t(n2) + t(n3)) / 2;
  Coeff = (V(n2) - V(n1)) / ((t(n2) - tMax)^2 - (t(n1) - tMax)^2);
else
  Val_1 = (V(n2) - V(n1)) / (V(n2) - V(n3));

  b = 2 * (t(n2) - t(n1) + Val_1 * (t(n3) - t(n2)));
  c = Val_1 * (t(n2)^2 - t(n3)^2) + t(n1)^2 - t(n2)^2;

  tMax = -c / b;
  Coeff = (V(n2) - V(n1)) / ((t(n2) - tMax)^2 - (t(n1) - tMax)^2);
  %arbitrary which formula to use:
  %Coeff = (V(n3) - V(n1)) / ((t(n3) - tMax)^2 - (t(n1) - tMax)^2);
end

MaxV = V(n2) - Coeff * (t(n2) - tMax)^2;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Freq = GetSpikeFrequency(Times, t)
if(length(t) == 0)
  Freq = 0;
  return
end
tHalf = .5 * (t(1) + t(end));
NumEvents = length(Times);
if(NumEvents == 0)
  Freq = 0;
  return
else
  %Check if there are no events in the second half of the experiment
  %  if so, presumably it just took a LONG time to settle down, so
  %  label the cell as NOT spiking
  if(length(find(Times > tHalf)) == 0)
    Freq = 0;
    return
  end
end

tElapsed = t(end) - t(1);
Freq = NumEvents / tElapsed * 1000;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutStruct = StructifyList(InList)
OutStruct.List = InList;
GoodInd = find(isfinite(InList));
if(length(GoodInd) > 1)
  InList = InList(GoodInd);
  OutStruct.Mean = mean(InList);
  OutStruct.StdDev = std(InList);
  OutStruct.Variance = OutStruct.StdDev^2;
  OutStruct.CoefOfVar = OutStruct.StdDev / OutStruct.Mean;
else
  OutStruct.Mean = 0;
  OutStruct.StdDev = 0;
  OutStruct.Variance = 0;
  OutStruct.CoefOfVar = 0;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [VMid, dVMid] = GetInteriorPoint(V, Deriv)
V = sort(V);
n = length(V);
while(V(n) > 1.1 * V(n-1))
  n = n - 1;
  if(n == 1)
    n = length(V);
    break
  end
end
VMid = V(n) / 3.0;

Deriv = sort(Deriv);
n = length(V);
while(Deriv(n) > 1.1 * Deriv(n-1))
  n = n - 1;
  if(n == 1)
    n = length(V);
    break
  end
end
dVMid = Deriv(n) / 3.0;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotGetSpikeTimes(t, V, Deriv, PlotSubject, LowCutoff, HighCutoff)
if(ischar(PlotSubject) & length(PlotSubject) > 0)
  TitleStr = ['V vs. dV/dt for ', PlotSubject];
else
  TitleStr = 'V vs. dV/dt';
end

[VMid, dVMid] = GetInteriorPoint(V, Deriv);
h = NamedFigure(TitleStr);
set(h, 'WindowStyle', 'docked');
hold off
plot(V, Deriv, 'b-');
hold on
plot(VMid, dVMid, 'ro', 'MarkerFaceColor', 'r');
xlabel('Voltage (mV)', 'FontSize', 18)
ylabel('dV/dt (mV/ms)', 'FontSize', 18)
title(RealUnderscores(TitleStr), 'FontSize', 18)
hold off

if(nargin < 6)
  return
end

if(ischar(PlotSubject) & length(PlotSubject) > 0)
  TitleStr = ['dV/dt for ', PlotSubject'];
else
  TitleStr = 'dV/dt';
end
h = NamedFigure(TitleStr);
set(h, 'WindowStyle', 'docked');
hold off
plot(t/1000, Deriv, 'b-');
hold on
plot([t(1), t(end)]/1000, [LowCutoff, LowCutoff], 'g-');
plot([t(1), t(end)]/1000, [HighCutoff, HighCutoff], 'g-');
xlabel('Time (s)', 'FontSize', 18)
ylabel('dV/dt (mV/ms)', 'FontSize', 18)
title(RealUnderscores(TitleStr), 'FontSize', 18)
hold off
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotGetSpikes(t, V, Spike, Burst, PlotSubject)
SpikeTimes = Spike.Times;
BurstTimes = Burst.Times;
BurstLen = Burst.Durations.List;

Top = max(V);
Bottom = min(V);
Delta = 0.25 * (Top - Bottom);
Bottom = Bottom - Delta;
Top = Top + Delta;

t = t / 1000;
SpikeTimes = SpikeTimes / 1000;
BurstTimes = BurstTimes / 1000;
BurstLen = BurstLen / 1000;

if(ischar(PlotSubject) & length(PlotSubject) > 0)
  TitleStr = ['Spike/Burst Detection for ', PlotSubject];
else
  TitleStr = 'Spike/Burst Detection';
end

h = NamedFigure(TitleStr);
set(h, 'WindowStyle', 'docked');
clf;

whitebg(h, 'k');
hold on;
%first draw red rectangle to signify burst times
NumBursts = length(BurstTimes);
for n = 1:NumBursts
  tLow = BurstTimes(n);
  tHigh = tLow + BurstLen(n);
  fill([tLow, tHigh, tHigh, tLow], [Bottom, Bottom, Top, Top], 'b');
end

%next overlay lines indicating spikes
NumSpikes = length(SpikeTimes);
for n=1:NumSpikes
  plot([SpikeTimes(n), SpikeTimes(n)], [Bottom,Top], 'r-');
  
  plot(Spike.MaxV.t(n)/1000, Spike.MaxV.V.List(n), 'go', 'MarkerSize', ...
       4, 'MarkerFaceColor', 'g');
  plot(Spike.PreMaxCurve.t(n)/1000, Spike.PreMaxCurve.V.List(n), 'go', ...
       'MarkerSize', 4, 'MarkerFaceColor', 'g');
  plot(Spike.PostMaxCurve.t(n)/1000, Spike.PostMaxCurve.V.List(n), 'go', ...
       'MarkerSize', 4, 'MarkerFaceColor', 'g');
  plot(Spike.PreMinV.t(n)/1000, Spike.PreMinV.V.List(n), 'co', ...
       'MarkerSize', 4, 'MarkerFaceColor', 'c');
  plot(Spike.PostMinV.t(n)/1000, Spike.PostMinV.V.List(n), 'co', ...
       'MarkerSize', 4, 'MarkerFaceColor', 'c');
end

%finally draw the voltage trace:
plot(t, V, 'w-');
xlabel('Time (s)', 'FontSize', 18)
ylabel('Voltage (mV)', 'FontSize', 18)
title(RealUnderscores(TitleStr), 'FontSize', 18)
hold off;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotVar = DoPlot(PlotSubject)
if(ischar(PlotSubject))
  PlotVar = true;
else
  PlotVar = PlotSubject;
end
return