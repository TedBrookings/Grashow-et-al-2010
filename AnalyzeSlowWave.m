function SlowWave = AnalyzeSlowWave(t, V, Spikes, PlotSubject)
if(nargin < 4)
  PlotSubject = false;
end
PlotPhases = true;
DoSpectrum = true;
DoScalogram = false;

t = t / 1000;  %convert to seconds
DeltaT = t(2) - t(1);
SampleLen = t(end) - t(1);

[VWave, BackInds, SpikeAmp] = RemoveSpikes(V, Spikes, DeltaT);

SlowWave = GetPhases(VWave, DeltaT, BackInds, SpikeAmp, ...
		     PlotPhases, PlotSubject);

SlowWave = CheckRealSlowWave(SlowWave, Spikes, VWave, t);

if(DoSpectrum)
  if(DoScalogram)
    [Spectrum.Freq, Spectrum.Amplitude, Spectrum.Phase] ...
	= Scalogram(V, t, PlotSubject);
    Spectrum.AvgPower = GetSpectrumPower(Spectrum, SlowWave, [], ...
					 PlotSubject, 'Average Power');
    SlowInds = GetSlowInds(Spikes, t);
    if(length(SlowInds) > 0)
      Spectrum.SlowPower = GetSpectrumPower(Spectrum, SlowWave, SlowInds, ...
					    PlotSubject, ...
					    'Slow-wave Power');
    end
  else
    %[Spectrum.Power, Spectrum.PowerConf, Spectrum.Freq] = ...
    %pmtm(V - mean(V), 1.25, [], 1.0 / DeltaT);
    nw = 0.5 * round(2.0 * SampleLen / 30.0);
    if(nw < 1)
      nw = 1;
    end
    NumTapers = floor(2 * nw - 1);
    MaxFreq = 30.0;
    if(size(V, 2) > 1)
      V = V';
    end
    [Spectrum.Freq, Spectrum.Power, DUMMY, DUMMY, DUMMY, ...
     f_res_diam, Spectrum.PowerConf] = ...
	coh_mt(DeltaT, V - mean(V), nw, NumTapers, MaxFreq, 0.68, 2);
    
    if(DoPlot(PlotSubject))
      PSlow = interp1(Spectrum.Freq, Spectrum.Power, SlowWave.Freq);
      TitleStr = 'Spectrogram';
      if(ischar(PlotSubject) & length(PlotSubject) > 0)
	TitleStr = [TitleStr, ' for ', PlotSubject];
      end
      h = NamedFigure(TitleStr);
      set(h, 'WindowStyle', 'docked');
      hold off
      %plot(Spectrum.Freq, sqrt(Spectrum.PowerConf(:,:,1)), 'r-');
      %hold on
      %plot(Spectrum.Freq, sqrt(Spectrum.PowerConf(:,:,2)), 'b-');
      plot(Spectrum.Freq, sqrt(Spectrum.Power), 'k-')
      hold on
      plot(SlowWave.Freq, sqrt(PSlow), 'go', 'MarkerFaceColor', 'g');
      xlim([0, 10])
      hold off
      title(RealUnderscores(TitleStr), 'FontSize', 18);
      xlabel('Frequency (Hz)', 'FontSize', 18);
      ylabel('Power', 'FontSize', 18);
    end
  end
  SlowWave.Spectrum = Spectrum;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [VWave, BackInds, SpikeAmp] = RemoveSpikes(V, Spikes, DeltaT)
VWave = V;
BackInds = [];
SpikeAmp = 0;
if(length(Spikes) == 0)  %Spikes weren't analyzed
  return
end
n1List = Spikes.n1List;
NumSpikes = length(n1List);
if(NumSpikes == 0)  %Spikes weren't found
  return
end
n2List = Spikes.n2List;
BackInds = zeros(size(n2List));
NumV = length(V);

DeltaList = [];
SkipLen = 2 * ceil(30.0e-3 / DeltaT);
n1Stop = SkipLen;
HalfSkipLen = SkipLen;  %  used to be: round(0.5 * SkipLen);
                       %  but this encourages expanding both ways
for SNum = 1:NumSpikes
  n1 = n1List(SNum);
  n2 = n2List(SNum);
  Delta = n2 - n1;

  MaxExpand = 3 * Delta;
  MaxBad = Delta;
  Expand = 0;
  BadMove = 0;
  DeltaV = abs(V(n1) - V(n2));
  SpikeAmp = SpikeAmp + DeltaV;
  BestDV = DeltaV;
  Bestn1 = n1;
  Bestn2 = n2;
  if(SNum == NumSpikes)
    n2Stop = NumV - SkipLen;
  else
    n2Stop = min([n1List(SNum + 1) - 1, NumV - SkipLen]);
  end
  while(DeltaV > 0.1 && Expand < MaxExpand)
    if(n1 <= n1Stop)
      ExpandType = 1;
      if(n2 >= n2Stop)
	break
      end
    elseif(n2 >= n2Stop)
      ExpandType = -1;
    else
      DVMinus = abs(V(n1-SkipLen) - V(n2));
      DVPlus = abs(V(n1) - V(n2+SkipLen));
      DVBoth = abs(V(n1-HalfSkipLen) - V(n2+HalfSkipLen));
      if(DVMinus < DVPlus)
	if(DVMinus < DVBoth)
	  ExpandType = -1;
	else
	  ExpandType = 0;
	end
      elseif(DVPlus < DVBoth)
	ExpandType = 1;
      else
	ExpandType = 0;
      end
    end
    
    switch(ExpandType)
     case -1, n1 = n1 - 1; Expand = Expand + 1;
     case 1, n2 = n2 + 1; Expand = Expand + 1;
     otherwise, n1 = n1 - 1; n2 = n2 + 1; Expand = Expand + 2;
    end

    DeltaV = abs(V(n1) - V(n2));
    if(DeltaV < BestDV)
      BestDV = DeltaV;
      Bestn1 = n1;
      Bestn2 = n2;
    else
      BadMove = BadMove + 1;
      if(BadMove > MaxBad)
	break
      end
    end
  end
  n1 = Bestn1;
  n2 = Bestn2;
  
  BackInds(SNum) = n2;
  
  n1Stop = max([n2 + 1, SkipLen]);
  Ind = n1:n2;
  VWave(Ind) = interp1([n1, n2], [V(n1), V(n2)], Ind);
end
SpikeAmp = SpikeAmp / NumSpikes;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SlowWave = GetPhases(VWave, DeltaT, BackInds, SpikeAmp, ...
			      PlotPhases, PlotSubject)

SigFact = 0.5;  %avg. amplitude must be larger than SigFact * SpikeAmp
SlowWave.Phases = [];
SlowWave.Amplitudes = [];
SlowWave.MinInds = [];

%First low-pass filter VWave (below 10 Hz) to eliminate noise
WFilt = (2 * DeltaT) * 10.0;
[B,A] = butter(1, WFilt, 'low');
VWave = filtfilt(B, A, VWave);

[SlowWave.Freq, SlowWave.Sigma, SlowWave.Corr] = ...
    GetSlowWave(DeltaT, VWave, PlotSubject);
Freq = SlowWave.Freq;
if(Freq <= 0)
  PlotWaveAnalysis(VWave, [], [], DeltaT, PlotSubject, PlotPhases);
  return
end

IndPeriod = round(1.0 / (DeltaT * Freq));
Phases = zeros(size(VWave));

NumWave = length(VWave);
LookFact = 0.1 + 0.4 * SlowWave.Corr;
LookAhead = round(LookFact * IndPeriod);

%Get a list of points smaller/larger than any point within distance LookAhead
[VWaveSort, SortInd] = sort(VWave);
MinInds = SortInd(1);
MaxInds = SortInd(end);
NumV = length(VWave);
m = NumV;
for n=2:NumV
  m = m - 1;
  Diff = min(abs(MinInds - SortInd(n)));
  if(Diff > LookAhead)
    MinInds = [MinInds, SortInd(n)];
  end
  Diff = min(abs(MaxInds - SortInd(m)));
  if(Diff > LookAhead)
    MaxInds = [MaxInds, SortInd(m)];
  end
end

%Remove minima that are too close to the back ends of spikes
if(length(BackInds) > 0)
  n = 1;
  while(n <= length(MinInds))
    m = MinInds(n);
    n2 = find(BackInds <= m, 1, 'last');
    if(length(n2) == 0)
      n = n + 1;
      continue
    end
    n2 = BackInds(n2);
    if(m - n2 <= 1)  %Remove this MinInd
      MinInds = [MinInds(1:(n-1)), MinInds((n+1):end)];
    else
      n = n + 1;
    end
  end
end

%Refine the list to only include local minima/maxima
n = 2;
while(n <= length(MinInds))
  Test = MinInds(n);
  I1 = max([1, Test - LookAhead]);
  I2 = min([NumV, Test + LookAhead]);
  LocalMin = min(VWave(I1:I2));
  if(LocalMin < VWave(Test))
    MinInds = [MinInds(1:(n-1)), MinInds((n+1):end)];
  else
    n = n + 1;
  end
end
MinInds = sort(MinInds);
n = 2;
while(n <= length(MaxInds))
  Test = MaxInds(n);
  I1 = max([1, Test - LookAhead]);
  I2 = min([NumV, Test + LookAhead]);
  LocalMax = max(VWave(I1:I2));
  if(LocalMax > VWave(Test))
    MaxInds = [MaxInds(1:(n-1)), MaxInds((n+1):end)];
  else
    n = n + 1;
  end
end
MaxInds = sort(MaxInds);
LocalMins = VWave(MinInds);
LocalMaxes = VWave(MaxInds);

if(length(MinInds) < 3 || length(MaxInds) < 3)
  PlotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, PlotSubject, PlotPhases);
  return
end

%now make sure they alternate
if(MinInds(1) < MaxInds(1))
  Current = -1;
else
  Current = 1;
end
nMin = 1;
nMax = 1;
Done = false;
while(~Done)
  if(Current > 0)
    if(nMax == length(MaxInds))
      Done = true;
      %prevent several mins in a row at the end:
      if(nMin < length(MinInds))
	[TempMin, TempInd] = min(VWave((MaxInds(end)+1):end));
	TempInd = TempInd + MaxInds(end);
	nMin = nMin - 1;
	MinInds = [MinInds(1:nMin), TempInd];
	LocalMins = [LocalMins(1:nMin), TempMin];
      end
      break;
    end
    if(MaxInds(nMax + 1) < MinInds(nMin))
      %Two maxes in a row.  Remove smaller
      if(LocalMaxes(nMax) > LocalMaxes(nMax + 1))
	nDel = nMax + 1;
      else
	nDel = nMax;
      end
      MaxInds = [MaxInds(1:(nDel - 1)), MaxInds((nDel + 1):end)];
      LocalMaxes = [LocalMaxes(1:(nDel - 1)), LocalMaxes((nDel + 1):end)];
    else
      nMax = nMax + 1;
      Current = -1;
    end
  else
    if(nMin == length(MinInds))
      Done = true;
      %prevent several maxes in a row at the end:
      if(nMax < length(MaxInds))
	[TempMax, TempInd] = max(VWave((MinInds(end)+1):end));
	TempInd = TempInd + MinInds(end);
	nMax = nMax - 1;
	MaxInds = [MaxInds(1:nMax), TempInd];
	LocalMaxes = [LocalMaxes(1:nMax), TempMax];
      end
      break;
    end
    if(MinInds(nMin + 1) < MaxInds(nMax))
      %Two mins in a row.  Remove bigger
      if(LocalMins(nMin) < LocalMins(nMin + 1))
	nDel = nMin + 1;
      else
	nDel = nMin;
      end
      MinInds = [MinInds(1:(nDel - 1)), MinInds((nDel + 1):end)];
      LocalMins = [LocalMins(1:(nDel - 1)), LocalMins((nDel + 1):end)];
    else
      nMin = nMin + 1;
      Current = 1;
    end
  end
end

if(length(MinInds) < 3 || length(MaxInds) < 3)
  PlotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, PlotSubject, PlotPhases);
  return
end

%ensure that each min is a local min between maxes, and vice versa
Changed = true;
while(Changed)
  Changed = false;
  if(MinInds(1) < MaxInds(1))
    Current = 1;
  else
    Current = -1;
  end
  nMin = 1;
  nMax = 1;
  Done = false;
  while(~Done)
    if(Current > 0)
      if(nMin == length(MinInds))
	Done = true;
	I2 = NumV;
      else
	I2 = MinInds(nMin + 1) - 1;
      end
      I1 = MinInds(nMin) + 1;
      [Val, Ind] = max(VWave(I1:I2));
      if(Val > LocalMaxes(nMax))
	%max out of place
	Changed = true;
	MaxInds(nMax) = Ind + I1 - 1;
	LocalMaxes(nMax) = Val;
      end
      nMin = nMin + 1;
      Current = -1;
    else
      if(nMax == length(MaxInds))
	Done = true;
	I2 = NumV;
      else
	I2 = MaxInds(nMax + 1) - 1;
      end
      I1 = MaxInds(nMax) + 1;
      [Val, Ind] = min(VWave(I1:I2));
      if(Val < LocalMins(nMin))
	%max out of place
	Changed = true;
	MinInds(nMin) = Ind + I1 - 1;
	LocalMins(nMin) = Val;
      end
      nMax = nMax + 1;
      Current = 1;
    end 
  end
end

if(length(MinInds) < 3 || length(MaxInds) < 3)
  PlotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, PlotSubject, PlotPhases);
  return
end

DebugRefine = false;
%Remove min/max pairs until the slow wave frequency seems to match
DeltaFreq = 1.0 / (DeltaT * (MinInds(end) - MinInds(1)));
Freq = (length(MinInds) - 1) * DeltaFreq;
DiffMaxes = GetDiffMeasure(LocalMaxes, 1, DebugRefine);
DiffMins = GetDiffMeasure(LocalMins, -1, DebugRefine);
if(DebugRefine)
  fprintf('CurrentFreq: %g, Actual: %g\n', Freq, SlowWave.Freq)
end

while(Freq - SlowWave.Freq > DeltaFreq || ...
      DiffMaxes >= 3 || DiffMins >= 3)
  if(DiffMaxes > DiffMins)
    [Dummy, RMaxInd] = min(LocalMaxes);
    IndMax = MaxInds(RMaxInd);
    Ind2 = find(MinInds > IndMax, 1);
    if(length(Ind2) == 0)
      RMinInd = length(LocalMins);
    else
      Ind1 = Ind2 - 1;
      if(Ind1 == 0 || LocalMins(Ind1) < LocalMins(Ind2))
	RMinInd = Ind2;
      else
	RMinInd = Ind1;
      end
    end
  else
    [Dummy, RMinInd] = max(LocalMins);
    IndMin = MinInds(RMinInd);
    Ind2 = find(MaxInds > IndMin, 1);
    if(length(Ind2) == 0)
      RMaxInd = length(LocalMaxes);
    else
      Ind1 = Ind2 - 1;
      if(Ind1 == 0 || LocalMaxes(Ind1) > LocalMaxes(Ind2))
	RMaxInd = Ind2;
      else
	RMaxInd = Ind1;
      end    
    end
  end
  
  MinInds = [MinInds(1:(RMinInd-1)), MinInds((RMinInd+1):end)];
  LocalMins = [LocalMins(1:(RMinInd-1)), LocalMins((RMinInd+1):end)];
  MaxInds = [MaxInds(1:(RMaxInd-1)), MaxInds((RMaxInd+1):end)];
  LocalMaxes = [LocalMaxes(1:(RMaxInd-1)), LocalMaxes((RMaxInd+1):end)];
  
  if(length(MinInds) < 3 || length(MaxInds) < 3)
    PlotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, PlotSubject, PlotPhases);
    return
  end

  DeltaFreq = 1.0 / (DeltaT * (MinInds(end) - MinInds(1)));
  Freq = (length(MinInds) - 1) * DeltaFreq;
  if(DebugRefine)
    fprintf('CurrentFreq: %g, Actual: %g\n', Freq, SlowWave.Freq)
  end
  DiffMaxes = GetDiffMeasure(LocalMaxes, 1, DebugRefine);
  DiffMins = GetDiffMeasure(LocalMins, -1, DebugRefine);
end

%MinInds and MaxInds are calculated
%  now calculate Phases:
if(MinInds(1) < MaxInds(1))
  Current = -1;
  CurrentMax = LocalMaxes(1);
  nMin = 1;
  nMax = 0;
else
  Current = 1;
  CurrentMin = LocalMins(1);
  nMin = 0;
  nMax = 1;
end
Done = false;
IndStart = 1;
Amplitudes = [];
while(~Done)
  if(Current > 0)  %increasing
    if(nMax > length(MaxInds))
      IndStop = NumV;
      Done = true;
    else
      CurrentMax = LocalMaxes(nMax);
      IndStop = MaxInds(nMax);
    end
    Ind = IndStart:IndStop;
    IndStart = IndStop + 1;
    Slope = 2.0 / (CurrentMax - CurrentMin);
    Offset = 1.0 + CurrentMin * Slope;
    Phases(Ind) = 0.5 * pi + real(asin(Slope * VWave(Ind) - Offset));
    nMin = nMin + 1;
    Current = -1;
  else   %decreasing
    if(nMin > length(MinInds))
      IndStop = NumV;
      Done = true;
    else
      CurrentMin = LocalMins(nMin);
      IndStop = MinInds(nMin);
    end
    Ind = IndStart:IndStop;
    IndStart = IndStop + 1;
    Slope = 2.0 / (CurrentMax - CurrentMin);
    Offset = 1.0 + CurrentMin * Slope;
    Phases(Ind) = 1.5 * pi - real(asin(Slope * VWave(Ind) - Offset));
    nMax = nMax + 1;
    Current = 1;
    if(nMax <= length(MaxInds))
      Amplitudes = [Amplitudes, 0.5 * (LocalMaxes(nMax) + CurrentMax) ...
		    - CurrentMin];
    end
  end
end

%fprintf('MeanAmp: %g, SpikeAmp: %g, CutOff:  %g\n', ...
%	mean(Amplitudes), SpikeAmp, SigFact * SpikeAmp)
if(mean(Amplitudes) < SigFact * SpikeAmp)
  PlotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, PlotSubject, PlotPhases);
  return
end
SlowWave.Phases = Phases;
SlowWave.Amplitudes = Amplitudes;
SlowWave.MinInds = MinInds;

PlotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, PlotSubject, PlotPhases);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, PlotSubject, ...
                          PlotPhases)
IntensePlot = false;
if(DoPlot(PlotSubject) & PlotPhases)
  TitleStr = 'Slow-wave Phase';
  if(ischar(PlotSubject) & length(PlotSubject) > 0)
    TitleStr = [TitleStr, ' for ', PlotSubject];
  end
  h = NamedFigure(TitleStr);
  set(h, 'WindowStyle', 'docked');
  clf
  hold on
  if(IntensePlot)
    t = 0;
    for n = 1:length(VWave)
      Blue = 0.5*(1.0 + cos(Phases(n)));
      Red = 1.0 - Blue;
      plot(t, VWave(n), '.', 'MarkerEdgeColor', [Red, 0, Blue], ...
	   'MarkerSize', 1)
      t = t + DeltaT;
    end
  else
    plot(DeltaT * (0:(length(VWave)-1)), VWave, 'b-');
  end
  if(length(MaxInds) >= 3 && length(MinInds) >= 3)
    plot(MaxInds * DeltaT, VWave(MaxInds), 'go');
    plot(MinInds * DeltaT, VWave(MinInds), 'gx');
  end
  hold off
  title(RealUnderscores(TitleStr), 'FontSize', 18);
  xlabel('Time (s)', 'FontSize', 18);
  ylabel('Voltage (mV)', 'FontSize', 18);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZWorst = GetWorstZ(X)
Z = zeros(size(X));
%fprintf('%s', 'Z:')
for n = 1:length(X)
  TempX = [X(1:(n-1)), X((n+1):end)];
  Z(n) = (X(n) - mean(TempX)) / std(TempX);
  %fprintf(' %g', Z(n))
end
%fprintf('\n')
%Z = zscore(X);

ZWorst = max(Z);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DiffMeasure = GetDiffMeasure(X, Type, DebugRefine)
%  Type is +1 for Maxes, -1 for Mins
d = pdist(X');
Z = linkage(d);
c = cluster(Z, 'maxclust', 2);

Ind1 = find(c == 1);
Mean1 = mean(X(Ind1));
Std1 = std(X(Ind1));
Ind2 = find(c == 2);
Mean2 = mean(X(Ind2));
Std2 = std(X(Ind2));

%look for abnormally small Maxes, abnormally large Mins
if(Mean1 * Type > Mean2 * Type)
  L1 = length(Ind1);
  L2 = length(Ind2);
else
  L1 = length(Ind2);
  L2 = length(Ind1);
end
if(L1 > L2 || L1 >= 3)  %could be abnormal in way we are concerned
  DiffScale = 2 * max(Std1, Std2);
  if(DiffScale < 1)
    DiffScale = 1;
  end
  DiffMeasure = abs(Mean2 - Mean1) / DiffScale;
else   %abnormally large Maxes or abnormally small Mins, don't care:
  DiffMeasure = 0;
end

if(DebugRefine)
  fprintf('M1: %g, S1: %g, M2: %g, S2: %g, Diff: %g\n', ...
          Mean1, Std1, Mean2, Std2, DiffMeasure)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MinInds, MaxInds] = CheckNew(VWave, MinInds, MaxInds, Type, ...
				       IndMin, IndMax)
if(nargin < 6)
  IndMin = length(MinInds);
  IndMax = length(MaxInds);
end
if(strcmp(Type, 'min'))
  if(IndMin < 2)
    return
  end
  Ind = MinInds(IndMin - 1):MinInds(IndMin);
  [Val, NewMaxInd] = max(VWave(Ind));
  NewMaxInd = NewMaxInd + Ind(1) - 1;
  if(MaxInds(IndMax) ~= NewMaxInd)
    MaxInds(IndMax) = NewMaxInd;
    [MinInds, MaxInds] = CheckNew(VWave, MinInds, MaxInds, 'max', ...
				  IndMin - 1, IndMax);
  end
else
  if(IndMax < 2)
    return
  end
  Ind = MaxInds(IndMax - 1):MaxInds(IndMax);
  [Val, NewMinInd] = min(VWave(Ind));
  NewMinInd = NewMinInd + Ind(1) - 1;
  if(MinInds(IndMin) ~= NewMinInd)
    MinInds(IndMin) = NewMinInd;
    [MinInds, MaxInds] = CheckNew(VWave, MinInds, MaxInds, 'min', ...
				  IndMin, IndMax - 1);
  end  
end

return
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SlowWave = CheckRealSlowWave(SlowWave, Spikes, VWave, t)
%First check to see if the slow-wave is significantly larger than
%  background noise

if(SlowWave.Freq <= 0)
  SlowWave.NumBurst = 0;
  SlowWave.NumOnlySpike = 0;
  SlowWave.NumNoSpike = 0;
  SlowWave.SpikesPerWave = [];
  return
end

%high-pass filter the waveform
SmoothTime = .1 * (1.0 / SlowWave.Freq);
n = ceil(SmoothTime / (t(2) - t(1)));
if(n < 16)
  n = 16;
end
[B,A] = butter(2, 2 / n, 'high');  %2/n because Nyquist rate is 1/2
VFilt = filtfilt(B, A, VWave);

%h = NamedFigure('poo');
%set(h, 'WindowStyle', 'docked');
%plot(t, VFilt)
VFilt = sort(abs(VFilt));
OneSigma = 0.682689492;
TwoSigma = 0.954499736;
Fuzz = 2 * VFilt(round(OneSigma * length(VFilt)));  %Factor of two
                                                    %makes two-sided

MedianAmp = median(SlowWave.Amplitudes);
SlowWave.Sigma = MedianAmp / Fuzz;
if(SlowWave.Sigma < 2)
  SlowWave.Amplitudes = [];
  SlowWave.MinInds = [];
  SlowWave.NumBurst = 0;
  SlowWave.NumOnlySpike = 0;
  SlowWave.NumNoSpike = 0;
  SlowWave.SpikesPerWave = [];
  return
end
  
%Now check to see if spikes are the cause of the slow-wave.
% If they are, the phase should be advanced (near trough of slow-wave)
% at the end of the spike.
if(length(Spikes) == 0)
  NumSpikes = 0;
else
  n1List = Spikes.n1List;
  NumSpikes = length(n1List);
end
NumWaves = length(SlowWave.MinInds) - 1;
if(NumSpikes == 0)
  SlowWave.NumBurst = 0;
  SlowWave.NumOnlySpike = 0;
  SlowWave.NumNoSpike = NumWaves;
  SlowWave.SpikesPerWave = zeros(NumWaves, 1);
  return
end
n2List = Spikes.n2List;
Phases = SlowWave.Phases;
MinInds = SlowWave.MinInds;

PhaseCutoff = 1.5 * pi;

NumBurst = 0;  %Increment if the phase is not advanced
NumOnlySpike = 0;  %Increment if the phase is advanced
NumNoSpike = 0;  %Increment if there is no spike at all
SpikesPerWave = zeros(NumWaves, 1);
for WaveNum = 1:NumWaves
  StartInd = MinInds(WaveNum);
  StopInd = MinInds(WaveNum + 1);
  SpikeInds = find(n1List > StartInd & n2List < StopInd);
  Num = length(SpikeInds);
  SpikesPerWave(WaveNum) = Num;
  if(Num == 0)
    NumNoSpike = NumNoSpike + 1;
  elseif(Num == 1)
    n2 = n2List(SpikeInds(end));
    Theta_n2 = Phases(n2);
    MinPhase = min(Phases(n2:StopInd-1));
    %fprintf('p_n2: %g, MinP: %g\n', Theta_n2, MinPhase))
    if(MinPhase > PhaseCutoff)
      NumOnlySpike = NumOnlySpike + 1;
    else
      NumBurst = NumBurst + 1;
    end
  else
    NumBurst = NumBurst + 1;
  end
end

SlowWave.NumBurst = NumBurst;
SlowWave.NumOnlySpike = NumOnlySpike;
SlowWave.NumNoSpike = NumNoSpike;
SlowWave.SpikesPerWave = SpikesPerWave;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Power = GetSpectrumPower(Spectrum, SlowWave, Indices, ...
				  PlotSubject, TitleStr)
if(length(Indices) == 0)
  Power = mean(Spectrum.Amplitude.^2, 2);
else
  Power = mean(Spectrum.Amplitude(:,Indices).^2, 2);
end

if(DoPlot(PlotSubject))
  PSlow = interp1(Spectrum.Freq, Power, SlowWave.Freq);
  if(ischar(PlotSubject) & length(PlotSubject) > 0)
    TitleStr = [TitleStr, ' for ', PlotSubject];
  end
  h = NamedFigure(TitleStr);
  set(h, 'WindowStyle', 'docked');
  hold off
  loglog(Spectrum.Freq, Power)
  hold on
  plot(SlowWave.Freq, PSlow, 'ro');
  hold off
  title(RealUnderscores(TitleStr), 'FontSize', 18);
  xlabel('Frequency (Hz)', 'FontSize', 18);
  ylabel('Power', 'FontSize', 18);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotVar = DoPlot(PlotSubject)
if(ischar(PlotSubject))
  PlotVar = true;
else
  PlotVar = PlotSubject;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SlowInds = GetSlowInds(Spikes, t)
n1List = Spikes.n1List;
n2List = Spikes.n2List;
NumT = length(t);
if(length(n1List) == 0)
  SlowInds = 1:NumT;
  return
end

if(n1List(1) > 1)
  SlowInds = 1:(n1List(1)-1);
else
  SlowInds = [];
end
SlowInds = [];

for n=2:length(n1List)
  StartInd = n2List(n-1) + 1;
  StopInd = n1List(n) - 1;
  Mid = round(0.5* (StartInd + StopInd));
  HalfRange = round(0.2 * 0.5 * (StopInd - StartInd));
  StartInd = Mid - HalfRange;
  StopInd = Mid + HalfRange;
  
  if(StopInd >= StartInd)
    SlowInds = [SlowInds, StartInd:StopInd];
  end
end

%if(n2List(end) < NumT)
%  SlowInds = [SlowInds, (n2List(end)+1):NumT];
%end

return