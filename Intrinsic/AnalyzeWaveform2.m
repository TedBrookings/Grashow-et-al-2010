function StructOut = AnalyzeWaveform2(t, v, varargin)
if(nargin < 3)
  PlotVar = false;
  NoShape = false;
  FirstOnly = false;
else
  PlotVar = varargin{1};
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

%May want to have different levels of analysis selectable

%T = (t(1):2:t(end));
%V = interp1(t, v, T);
%whos

%clear T V;
%StructOut = AnalyzeWaveform(t, V, LowThresh, HighThresh, PlotVar)
%analyzes a single voltage waveform, looking for spikes
%    and bursts, and calculating relevant frequencies.
%
%  INPUT PARAMETERS:
%   -t is time in ms
%   -V is voltage in mV
%    OPTIONAL:
%     -LowThresh sets the lower threshold (in dV/dt) for spikes
%     (defaults to 0)
%     -HighThresh sets the upper threshold.  If both LowThresh and
%       HighThresh are zero, they are determined automatically
%     -PlotVar should be set to 1 to produce plots of process
%  OUTPUT PARAMETERS:
%   -StructOut.Spike:  structure with spike information
%      -Spike.Freq is overall spiking frequency (in Hz)
%      -Spike.Times is list of spike times (in ms)
%      -Spike.Intervals is a list of interspike intervals (in ms)
%      -Spike.Frequencies is a list of instantaneous frequencies (in Hz)
%      -Spike.Shape is a structure with spike-shape information
%         -Shape.V is a list of voltages describing the shape (in mV)
%         -Shape.t is a corresponding list of times (in ms)
%           relative to time of maximum voltage
%         -Shape.VMaxes is a list of maximum voltages (in mV)
%         -Shape.VHalfs is a list of voltages between peak and trough (in mV)
%         -Shape.THalf1 is a list of times where voltage crosses
%           VHalfs before the maximum (in mV)
%         -Shape.THalf1 is a list of times where voltage crosses
%           VHalfs after the maximum (in mV)
%         -Shape.PreMins is a list of typical voltages before a spike (in mV)
%         -Shape.PostMins is a list of typical voltages after a spike (in mV)
%   -StructOut.Burst:  structure with burst information
%      -Burst.Freq is burst frequency (in Hz)
%      -Burst.SpikeFreq is within-burst spike frequency (in Hz)
%      -Burst.DutyCycle is the average burst duration/period
%      -Burst.Times is a list of burst times (in ms)
%      -Burst.Durations is a list of burst durations (in ms)
%      -Burst.NumSpikes is a list of spikes per burst
%      -Burst.SpikeFrequencies is a list of spike frequencies (in Hz)
%      -Burst.InterBurstIntervals is a list of inter-burst
%       intervals (in ms)
%   -StructOut.MedianV:  the median of the voltage trace.  If the
%      cell is silent, it should be the resting potential,
%      otherwise, who knows...
%If a feature is not detected, relevant frequencies are set to
%  zero, and relevant lists are empty

%First get the spike times
Spike = GetSpikeTimes(t, v, PlotVar, NoShape, FirstOnly);

%Next get the overall spike frequency
Spike.Freq = GetSpikeFrequency(Spike.Times, t);

%Now look for bursts, as deviations from overall spike frequency
GapFact = 2;
GapAvgFact = 2 * GapFact / (GapFact + 1);
Burst = GetBurstTimes(Spike, Spike.Freq, t, GapAvgFact);
%Find burst frequency, and within-burst spike frequency
Burst = GetBurstFreq(Burst);

%repeat to refine burst detection
NumRefine = 2;
for n = 1:NumRefine
  if(length(Burst.Times) == 0)
    break;
  end
  %Now look for bursts by finding deviations from within-burst frequency
  Burst = GetBurstTimes(Spike, Burst.SpikeFreq, t, GapFact);
  %Find updated burst-frequency and within-burst spike frequency
  Burst = GetBurstFreq(Burst);
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
Burst.NumSpikes = StructifyList(Burst.NumSpikes);
Burst.InterBurstIntervals = StructifyList(Burst.InterBurstIntervals);
Burst.SpikeFrequencies = StructifyList(Burst.SpikeFrequencies);

StructOut.Spike = Spike;
StructOut.Burst = Burst;
StructOut.MedianV = median(v);

if(PlotVar)
  PlotGetSpikes(t, v, Spike, Burst);
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Spike = GetSpikeTimes2(t, V, varargin)
if(nargin < 3)
  PlotVar = false;
  NoShape = false;
else
  PlotVar = varargin{1};
  if(nargin < 4)
    NoShape = false;
  else
    NoShape = varargin{2};
  end
end

DeltaT = t(2) - t(1);

Deriv = diff(V) / DeltaT;
if(NoShape)
  K = [];
else
  K = diff(Deriv) / DeltaT;
end
Deriv = .5 * (Deriv(2:end) + Deriv(1:(end-1)));
if(~NoShape)
  K = abs(K) .* (1 + Deriv.^2).^-1.5;
end

V = V(2:(end-1));
t = t(2:(end-1));

NumV = length(V);

VArr = [median(V), max(V)];
VMid = .5 * VArr(1) + .5 * VArr(2);

n1List = [];
n2List = [];

[Num, x] = hist(V, 100);
[Max, MaxInd] = max(Num);
VBot = x(MaxInd);


%VBot = x(HistInd(end) + 1)
VTop = max(V);
VMid = .5 * (VBot + VTop);
%VRange = VTop - VBot;
%dVRange = max(Deriv) - min(Deriv);


n = 1;
while(n <= NumV)
  if(V(n) < VMid)
    n = n + 1;
  else
    if(Deriv(n) < 0)
      n1List = [];
      n2List = [];
      break
    end
    n1 = n - 1;
    while(n1 > 0)
      if(Deriv(n1) < 0)
	break
      else
	n1 = n1 - 1;
      end
    end
    n1 = n1 + 1;
    n1List = [n1List, n1];
    
    n2 = n + 1;
    StopBelowVMid = true;
    if(StopBelowVMid)
      while(n2 <= NumV)
	if(V(n2) > VMid)
	  n2 = n2 + 1;
	else
	  break
	end
      end
    else
      while(n2 <= NumV)
	if(Deriv(n2) >= 0)
	  n2 = n2 + 1;
	else
	  break
	end
      end
    end
    while(n2 <= NumV)
      if(Deriv(n2) < 0)
	n2 = n2 + 1;
      else
	break
      end
    end
    n = n2;
    n2 = n2 - 1;
    n2List = [n2List, n2];
  end
end

%now spikes are all bracketed between n1List and n2List

Spike = GetSpikeShape(n1List, n2List, t, V, Deriv, K, NoShape);

if(PlotVar)
  PlotGetSpikeTimes(t, V, Deriv);
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
function Spike = GetSpikeTimes(t, V, varargin)
if(nargin < 3)
  PlotVar = false;
  NoShape = false;
  FirstOnly = false;
else
  PlotVar = varargin{1};
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

Deriv = diff(V) / DeltaT;
if(size(Deriv, 1) > 1)
  Deriv = [0; .5 * (Deriv(2:end) + Deriv(1:(end-1))); 0];
else
  Deriv = [0, .5 * (Deriv(2:end) + Deriv(1:(end-1))), 0];
end

AutoCutoffs = true;
if(AutoCutoffs)
  %Get typical "noise" DerivLevels
  [Num, x] = hist(V, 100);
  %HistInd = find(Num > mean(Num));
  [MaxNum, MaxInd] = max(Num);
  HistInd = find(V <= x(MaxInd + 1));
  [max(V(HistInd)), min(V(HistInd))];
  [DerivLow, DerivHigh] = GetDVRange(Deriv, HistInd);
  LowCutoff = DerivLow * 3.5;
  HighCutoff = DerivHigh * 3.5;
else
  LowCutoff = -2;
  HighCutoff = 2;
end


n1List = [];
n2List = [];
n = 2;
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
    while(n2 <= NumV & t(n2) - t(n) < MaxTimeWidth)
      if(Deriv(n2) > LowCutoff)
	if(Deriv(n2) > HighCutoff)  %If we've gone back above cutoff,
	  n = n2;                   %reset n
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
      if(Deriv(n1) < 0)
	break
      else
	n1 = n1 - 1;
      end
    end
    n1 = n1 + 1;
    n1List = [n1List, n1];
    
    while(n2 <= NumV)
      if(Deriv(n2) > 0)
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

Spike = GetSpikeShape(n1List, n2List, t, V, Deriv, NoShape);

if(PlotVar)
  PlotGetSpikeTimes(t, V, Deriv, LowCutoff, HighCutoff);
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
function varargout = GetDVRange(Deriv, Ind)
DVals = sort(Deriv(Ind));
Num = length(DVals);
LowBall = DVals(round(.025 * Num));
HighBall = DVals(round(.975 * Num));
if(nargout == 1)
  DVRange = HighBall - LowBall;
  varargout = {DVRange};
else
  varargout = {LowBall, HighBall};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Spike = GetSpikeShape(n1List, n2List, t, V, Deriv, NoShape)
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
  
  SmoothTime = 2;
  DeltaT = t(2) - t(1);
  if(DeltaT < SmoothTime)
    n = 2;
    while(t(n) - t(1) < SmoothTime)
      n = n + 2;
    end
    [B,A] = butter(4, 2 / n,'low');
    VSmooth = filtfilt(B, A, V);
  else
    VSmooth = t;
  end
  DS = diff(VSmooth) / DeltaT;
  DS = .5 * (DS(2:end) + DS(1:(end-1)));
  K = diff(DS) / DeltaT;
  DS = .5 * (DS(2:end) + DS(1:(end-1)));
  K = K .* (1 + DS.^2).^-1.5;
  if(size(K,1) > 1)
    K = [0; 0; K; 0; 0];
  else
    K = [0, 0, K, 0, 0];
  end
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
  [PreMaxK, tPreMaxK] = GetExtremum(K, t, n1, nMaxDV, 'max', true);
  VPreMaxK = interp1(t, V, tPreMaxK);
  [PostMaxK, tPostMaxK] = GetExtremum(K, t, nMinDV, n2, 'max', true);
  VPostMaxK = interp1(t, V, tPostMaxK);
  
  %Find minimum voltage before and after spike
  if(n1 > 1)
    n1 = n1 - 1;
  end
  if(n2 < length(V))
    n2 = n2 + 1;
  end
  [PreMinV, tPreMin] = GetExtremum(V, t, n1, n1+3, 'min');
  [PostMinV, tPostMin] = GetExtremum(V, t, n2-3, n2, 'min');
  
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
function Burst = GetBurstTimes(Spike, SpikeFreq, t, GapFact)
if(nargin < 4)
  GapFact = 4/3;
elseif(GapFact == 0)
  GapFact = 4/3;
end

Burst.Times = [];
Burst.Durations = [];
Burst.NumSpikes = [];
Burst.InterBurstIntervals = [];
Burst.DutyCycle = 0;

NumSpikes = length(Spike.Times);
if(NumSpikes == 0)
  return
end
%To be the beginning of a burst, require a leading time gap between
%  spikes that is GapFact times the inter-spike intervals
GapCutoff = GapFact * 1000 / SpikeFreq;
%To be a burst, there must be at least two spikes
SpikeCutoff = 2;


DeltaT = [Spike.Times(1) - t(1), Spike.Intervals];
BurstInd = find(DeltaT > GapCutoff);
if(length(BurstInd) == 0)
  return
end

%The length of the burst (in ms) is the time from the first spike
%   in the burst to the last spike in the burst:

BurstSpikes = NumSpikes - BurstInd(end);
BurstLen = Spike.Times(end) - Spike.Times(BurstInd(end));
if(length(BurstInd) > 1)
  BurstSpikes = [BurstInd(2:end) - BurstInd(1:(end-1)), BurstSpikes];
  BurstLen = [Spike.Times(BurstInd(2:end) - 1) ...
	      - Spike.Times(BurstInd(1:(end-1))), BurstLen];
end


Ind = find(BurstSpikes >= SpikeCutoff);
BurstInd = BurstInd(Ind);
BurstLen = BurstLen(Ind);
BurstSpikes = BurstSpikes(Ind);

if(length(BurstInd) == 0)
  return
end

MeanSpikes = sum(BurstSpikes) / length(BurstSpikes);
SigmaSpikes = (sum((BurstSpikes - MeanSpikes).^2)/ length(BurstSpikes))^.5;
if(SigmaSpikes > .2 * MeanSpikes)  %really just some Poisson train,
                                   %no actual bursts
  return
end

Burst.Times = Spike.Times(BurstInd);
Burst.Durations = BurstLen;
Burst.NumSpikes = BurstSpikes;

Burst.InterBurstIntervals = Burst.Times(2:end) - Burst.Times(1:(end-1)) ...
    - Burst.Durations(1:(end-1));

Burst.DutyCycle = sum(Burst.Durations(1:(end-1))) ...
    / (Burst.Times(end) - Burst.Times(1));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Burst = GetBurstFreq(Burst)
if(length(Burst.Times) < 2)
  Burst.Freq = 0;
  Burst.SpikeFreq = 0;
  Burst.SpikeFrequencies = [];
  return
end

DeltaT = Burst.Times(2:end) - Burst.Times(1:(end-1));
Burst.Freq = 1000 * length(DeltaT) / sum(DeltaT);

Weight = (Burst.NumSpikes - 1) / sum(Burst.NumSpikes - 1);
Burst.SpikeFrequencies = 1000 * (Burst.NumSpikes - 1) ./ Burst.Durations;
Burst.SpikeFreq = sum(Weight .* Burst.SpikeFrequencies);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutStruct = StructifyList(InList)
OutStruct.List = InList;
if(length(InList) > 1)
  OutStruct.Mean = mean(InList);
  OutStruct.Variance = sum((InList - OutStruct.Mean).^2) ...
      / (length(InList) - 1);
  OutStruct.StdDev = OutStruct.Variance^.5;
  OutStruct.CoefOfVar = OutStruct.StdDev / OutStruct.Mean;
else
  OutStruct.Mean = 0;
  OutStruct.Variance = 0;
  OutStruct.StdDev = 0;
  OutStruct.CoefOfVar = 0;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotGetSpikeTimes(t, V, Deriv, LowCutoff, HighCutoff)
h = figure;
set(h, 'WindowStyle', 'docked');
plot(V, Deriv, 'b-');

VArr = [median(V), max(V)];
VMid = .5 * VArr(1) + .5 * VArr(2);
hold on
plot(VMid, 0, 'ro');
hold off

h = figure;
set(h, 'WindowStyle', 'docked');
plot(t, Deriv, 'b-');
if(nargin == 5)
  hold on
  plot([t(1), t(end)], [LowCutoff, LowCutoff], 'g-');
  plot([t(1), t(end)], [HighCutoff, HighCutoff], 'g-');
  hold off
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotGetSpikes(t, V, Spike, Burst)
SpikeTimes = Spike.Times;
BurstTimes = Burst.Times;
BurstLen = Burst.Durations.List;

Top = max(V);
Bottom = min(V);
Delta = Top - Bottom;
Bottom = Bottom - 1.1 * Delta;
Top = Top + 1.1 * Delta;

%t = t / 1000;
%SpikeTimes = SpikeTimes / 1000;
%BurstTimes = BurstTimes / 1000;
%BurstLen = BurstLen / 1000;

h = figure;
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
  plot(Spike.MaxV.t(n), Spike.MaxV.V.List(n), 'g.');
  plot(Spike.PreMaxCurve.t(n), Spike.PreMaxCurve.V.List(n), 'g.');
  plot(Spike.PostMaxCurve.t(n), Spike.PostMaxCurve.V.List(n), 'g.');
end

%finally draw the voltage trace:
plot(t, V, 'w-');

hold off;
title('Spikes and Bursts')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotGetFreq(DeltaT, DenseLen, Densities)
DeltaTAvg = zeros(size(Densities));
for n = 1:length(Densities)
  DeltaTAvg(n) = sum(DeltaT(n:(n+DenseLen - 1))) / DenseLen;
end

h = figure;
set(h, 'WindowStyle', 'docked');
clf;
plot(DeltaTAvg, Densities);
return