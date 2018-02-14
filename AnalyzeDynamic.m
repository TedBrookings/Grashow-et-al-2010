function StructOut = AnalyzeDynamic(t, V0, V1, StabilizationTime, varargin)
% StructOut = AnalyzeDynamic(t, V0, V1, StabilizationTime, ...
%                            LowThresh, HighThresh, PlotSubject)
% examines two voltage traces (presumably from a dynamic clamp experiment)
% and calculates various properties of each trace and their mutual interactions
%
%  INPUT PARAMETERS:
%   -t is time in ms
%   -V0 is voltage of first cell (in mV)
%   -V1 is voltage of second cell
%   -StabilizationTime is duration at start of data that should be IGNORED
%    OPTIONAL:
%     -LowThresh sets the lower threshold (in dV/dt) for spikes
%     (defaults to 0)
%     -HighThresh sets the upper threshold
%     -PlotSubject should be set to true[false] to produce[suppress]
%       plots of waveforms/analysis.  Alternatively, it can be set
%       to a string to aid it titling plots (e.g. 'Exp #71')
%       PlotVar defaults to false
%  OUTPUT PARAMETERS:
%    -StructOut.Cell0:  structure with information about cell 0
%       returned from AnalyzeWaveform.m.  It should be
%       straightforward, but type 'help AnalyzeWaveform' for info.
%    -StructOut.Cell1:  structure with information about cell 1
%    -StructOut.SlowWave:  structure with information about
%     slow-wave behaviour of the combined system
%       -SlowWave.Freq is the frequency of the dominant slow-wave
%        component (in Hz)
%       -SlowWave.Sigma is a (very crude) measure of the importance
%        of the slow-wave frequency in the power spectrum
%    -StructOut.HalfCenter: structure with information about
%     half-center oscillations
%       -HalfCenter.Freq is the half-center frequency (in Hz)
%       -HalfCenter.Phase is the phase that cell 1 lags behind cell0
%       -HalfCenter.PhaseSigma is the standard deviation of the phases
%       -HalfCenter.ExclusionFact is the tendancy of bursts to
%        avoid each other.  0 is random, 1 is maximum, and -1 is minimum.
%       -HalfCenter.Periods is a list of full half-center periods (in ms)
%
% If a feature is not detected, relevant frequencies are set to
% zero, and relevant lists are empty

%close all;
switch nargin
 case 5,
  LowThresh = varargin{1};
  HighThresh = 0;
  PlotSubject = false;
 case 6,
  LowThresh = varargin{1};
  HighThresh = varargin{2};
  PlotSubject = false;
 case 7,
  LowThresh = varargin{1};
  HighThresh = varargin{2};
  PlotSubject = varargin{3};
 otherwise
  error(sprintf(['Incorrect number of arguments to AnalyzeDynamic.\n', ...
		 'Type "help AnalyzeDynamic" for details.']));
end

n1 = find(t > t(1) + StabilizationTime * 1000, 1);
t = t(n1:end);
V0 = V0(n1:end);
V1 = V1(n1:end);

if(DoPlot(PlotSubject))
  if(ischar(PlotSubject))
    PlotS0 = [PlotSubject, ' Cell0'];
    PlotS1 = [PlotSubject, ' Cell1'];
  else
    PlotS0 = 'Cell0';
    PlotS1 = 'Cell1';
  end
else
  PlotS0 = false;
  PlotS1 = false;
end
Cell0 = AnalyzeWaveform3(t, V0, PlotS0);
Cell1 = AnalyzeWaveform3(t, V1, PlotS1);

DeltaT = (t(2) - t(1)) / 1000;
[SlowWave.Freq, SlowWave.Sigma] = GetSlowWave(DeltaT, V0, V1, PlotSubject);
SpikeFreqAvg = GetSpikeFreqAvg(Cell0.Spike.Freq, Cell1.Spike.Freq);

HalfCenter = CheckHalfCenter(Cell0, Cell1, PlotSubject);
HalfCenter.Periods = StructifyList(HalfCenter.Periods);

disp(sprintf('SpikeFreq0 = %gHz, BurstFreq0 = %gHz', ...
	     Cell0.Spike.Freq, Cell0.Burst.Freq))
disp(sprintf('SpikeFreq1 = %gHz, BurstFreq1 = %gHz', ...
	     Cell1.Spike.Freq, Cell1.Burst.Freq))
disp(sprintf('SlowFreq = %gHz, NumSigma = %g', ...
	     SlowWave.Freq, SlowWave.Sigma))
disp(sprintf('Half-CenterFreq = %gHz, Exclusion Factor = %g', ...
	     HalfCenter.Freq, HalfCenter.ExclusionFact))
disp(sprintf('Half-CenterPhase = %g +- %g', ...
	     HalfCenter.Phase, HalfCenter.PhaseSigma))

StructOut.Cell0 = Cell0;
StructOut.Cell1 = Cell1;
StructOut.SlowWave = SlowWave;
StructOut.HalfCenter = HalfCenter;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SpikeFreqAvg = GetSpikeFreqAvg(SpikeFreq0, SpikeFreq1)
if(SpikeFreq0 == 0)
  if(SpikeFreq1 == 0)
    SpikeFreqAvg = 100;
  else
    SpikeFreqAvg = SpikeFreq1;
  end
else
  if(SpikeFreq1 == 0)
    SpikeFreqAvg = SpikeFreq0;
  else
    SpikeFreqAvg = 0.5 * (SpikeFreq0 + SpikeFreq1);
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HalfCenter = CheckHalfCenter(Cell0, Cell1, PlotSubject)
if(length(Cell0.Burst.Times) == 0)
  Cell0 = FakeBurst(Cell0);
end
if(length(Cell1.Burst.Times) == 0)
  Cell1 = FakeBurst(Cell1);
end

if(length(Cell0.Burst.Times) == 0 | length(Cell1.Burst.Times) == 0)
  HalfCenter.Freq = 0;
  HalfCenter.Phase = 0;
  HalfCenter.PhaseSigma = Inf;
  HalfCenter.ExclusionFact = 0;
  HalfCenter.Periods = [];
  return
end

%We want a train of alternating bursts (perhaps the cycle has two bursts
%  for cell 0, followed by one for cell 1, etc.  We want the first in
%  each cycle).
MixedTrain = GetBurstTrain(Cell0.Burst.Times, Cell1.Burst.Times);
NumTrain = length(MixedTrain);
if(NumTrain < 3)
  HalfCenter.Freq = 0;
  HalfCenter.Phase = 0;
  HalfCenter.PhaseSigma = Inf;
  HalfCenter.ExclusionFact = 0;
  HalfCenter.Periods = [];
  return
end

Intervals = MixedTrain(3:end) - MixedTrain(1:(end-2));
Phases = 2 * pi * (MixedTrain(2:(end - 1)) - MixedTrain(1:(end-2))) ...
			      ./ Intervals;
%Make all phases relative to the same starting cell: 
Phases(2:2:end) = 2 * pi - Phases(2:2:end);

NumIntervals = length(Intervals);

HalfCenter.Freq = 1000 * NumIntervals / sum(Intervals);
PhaseSum = sum(exp(i * Phases)) / NumIntervals;
HalfCenter.Phase = angle(PhaseSum);
HalfCenter.PhaseSigma = (-2 * log(abs(PhaseSum))).^.5;

HalfCenter.ExclusionFact = GetExclusionFact(Cell0.Burst.Times, ...
					    Cell0.Burst.Durations.List, ...
					    Cell1.Burst.Times, ...
					    Cell1.Burst.Durations.List);

HalfCenter.Periods = Intervals;

if(DoPlot(PlotSubject))
  PlotCheckHalfCenter(Intervals, Phases, PlotSubject);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cell = FakeBurst(Cell)
if(~(Cell.Spike.Freq > 0 & length(Cell.Spike.Times) >= 4))
  return
end
Cell.Burst.Times = Cell.Spike.Times(1:(end-1));

%Pick the BurstLen to be 1/4 the median inter-spike interval
DeltaT = sort(Cell.Spike.Times(2:end) - Cell.Burst.Times);
Duration = .25 * DeltaT(round(.5 * end));
Cell.Burst.Durations.List = repmat(Duration, size(Cell.Burst.Times));
Cell.Burst.NumSpikes.List = repmat(1, size(Cell.Burst.Times));

%Make the burst intervals symmetric about the spikes
Cell.Burst.Times = Cell.Burst.Times - Cell.Burst.Durations.List / 2;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MixedTrain = GetBurstTrain(BurstTimes0, BurstTimes1)
MixedTrain = [];
n0 = 1;
n1 = 1;

Len0 = length(BurstTimes0);
Len1 = length(BurstTimes1);

while(BurstTimes0(n0) == BurstTimes1(n1))
  n0 = n0 + 1;
  n1 = n1 + 1;
  if(n0 > Len0 | n1 > Len1)
    return
  end
end
if(BurstTimes0(1) < BurstTimes1(1))
  MixedTrain = [BurstTimes0(1), BurstTimes1(1)];
  Last = 1;
else
  MixedTrain = [BurstTimes1(1), BurstTimes0(1)];
  Last = 0;
end
n0 = n0 + 1;
n1 = n1 + 1;

while(1)
  if(Last == 0)
    if(n1 > Len1)
      return
    end
    while(BurstTimes1(n1) <= MixedTrain(end))
      n1 = n1 + 1;
      if(n1 > Len1)
	return
      end
    end
    MixedTrain = [MixedTrain, BurstTimes1(n1)];
    n1 = n1 + 1;
    Last = 1;
  else
    if(n0 > Len0)
      return
    end
    while(BurstTimes0(n0) <= MixedTrain(end))
      n0 = n0 + 1;
      if(n0 > Len0)
	return
      end
    end
    MixedTrain = [MixedTrain, BurstTimes0(n0)];
    n0 = n0 + 1;
    Last = 0;
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExclusionFact = GetExclusionFact(BurstTimes0, BurstLen0, ...
					  BurstTimes1, BurstLen1)
End0 = BurstTimes0 + BurstLen0;
End1 = BurstTimes1 + BurstLen1;
Ind0 = 1;
Ind1 = 1;
Num0 = length(BurstTimes0);
Num1 = length(BurstTimes1);
Overlap = 0;

while(Ind0 <= Num0 & Ind1 <= Num1)
  NewOverlap = min([End0(Ind0), End1(Ind1)]) ...
      - max([BurstTimes0(Ind0), BurstTimes1(Ind1)]);
  if(NewOverlap > 0)
    Overlap = Overlap + NewOverlap;
  end
  if(End0(Ind0) < End1(Ind1))
    Ind0 = Ind0 + 1;
  else
    Ind1 = Ind1 + 1;
  end
end

TotalTime = max([End0(end), End1(end)]) ...
    - min([BurstTimes0(1), BurstTimes1(1)]);

Len0 = sum(BurstLen0);
Len1 = sum(BurstLen1);
if(Len0 > Len1)
  DeltaLen = TotalTime - Len0;
  MinLen = Len1;
else
  DeltaLen = TotalTime - Len1;
  MinLen = Len0;
end

if(DeltaLen > MinLen)
  AvgOverlap = .5 * MinLen^2 / DeltaLen;
else
  AvgOverlap = MinLen - .5 * DeltaLen;
end

MinOverlap = Len0 + Len1 - TotalTime;
if(MinOverlap < 0)
  MinOverlap = 0;
end
MaxOverlap = MinLen;

if(Overlap < AvgOverlap)
  ExclusionFact = (AvgOverlap - Overlap) / (AvgOverlap - MinOverlap);
else
  ExclusionFact = (AvgOverlap - Overlap) / (MaxOverlap - AvgOverlap);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NumSig = GetPowerSigma(f, SigPower, Freq)
if(Freq == 0)
  NumSig = 0;
  return
end

FreqPow = interp1(f, SigPower, Freq);

SigPower = SigPower(find(SigPower ~= 0));
SigPower = sort(SigPower);
SigPower = SigPower(1:ceil(length(SigPower) * .9));
NumInd = length(SigPower);
PowBase = sum(SigPower) / NumInd;
PowSigma = sum( (SigPower - PowBase).^2 / NumInd ).^.5;

NumSig = (FreqPow - PowBase) / PowSigma;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutStruct = StructifyList(InList)
OutStruct.List = InList;
if(length(InList) > 1)
  OutStruct.Mean = sum(InList) / length(InList);
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
function PlotCheckHalfCenter(Intervals, Phases, PlotSubject)
Intervals = Intervals / 1000;  %convert from ms to s
if(ischar(PlotSubject) & length(PlotSubject) > 0)
  TitleStr = ['Burst Intervals for ', PlotSubject];
else
  TitleStr = 'Burst Intervals';
end
h = NamedFigure(TitleStr);
set(h, 'WindowStyle', 'docked');
hold off;
NumBars = length(Intervals) / 20;
if(NumBars < 20)
  NumBars = 20;
end
hist(Intervals, NumBars);
xlabel('Interval (s)', 'FontSize', 18);
ylabel('Number', 'FontSize', 18);
title(RealUnderscores(TitleStr), 'FontSize', 18);

if(ischar(PlotSubject) & length(PlotSubject) > 0)
  TitleStr = ['Burst Phases for ', PlotSubject];
else
  TitleStr = 'Burst Phases';
end
h = NamedFigure(TitleStr);
set(h, 'WindowStyle', 'docked');
hold off;
hist(Phases, NumBars);
xlim([0 2*pi])
xlabel('Phase', 'FontSize', 18);
ylabel('Number', 'FontSize', 18);
title(RealUnderscores(TitleStr), 'FontSize', 18);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotVar = DoPlot(PlotSubject)
if(ischar(PlotSubject))
  PlotVar = true;
else
  PlotVar = PlotSubject;
end
return