function [T, V, varargout] = CleanDCC(t, v, varargin)
% [T, V, Arr1, Arr2, ..., DCC_Info] = CleanDCC(t, v, arr1, arr2, ...)
% Resamples a voltage trace to a courser sample, in order to reduce
% noise introduced by DCC.
% DCC_Freq is calculated by looking for the dominant fast frequency
% of abs(dv/dt)  (because voltage shouldn't change much during the
% DCC period).
%  INPUTS:
%   -t: Array of times
%   -v: Num_t x NumTraces array of voltages
%   -arr1: (OPTIONAL) Num_t x NumTraces array of voltages/currents.
%  OUTPUTS:
%   -T: Array of times
%   -V: Num_T x NumTraces array of voltages
%   -Arr1: (OPTIONAL) Num_T x NumTraces array of currents, only
%       passed back if arr1 is passed in as non-empty array.
%   -DCC_Info: (OPTIONAL) A structure containing DCC information
%      .DCC_Freq:  The frequency of DCC sampling (in kHz)
%      .DCC_Power: The spectrum power at the DCC_Freq
%      .f:         An array of spectrum frequencies
%      .Power:     An Array of spectrum powers

DebugClean = false;

if(nargin < 2)
  error('Incorrect number of input arguments')
end
DeltaArgs = nargout - nargin;
if(DeltaArgs < 0 || DeltaArgs > 1)
  error('Invalid combination of in/out arguments.  Run "help CleanDCC"')
end

%First calculate the derivative and fast-frequency correlogram
NumAutoCorr = find(t - t(1) > 50, 1);  %Specify length
[AutoCorr, DV] = GetAutoCorr(t, v, NumAutoCorr);

%Next find the first maximum of autocorrelation (First guess for
%  DCC period)
[DCC_Ind, DCC_Corr] = GetAutoCorrMax(AutoCorr, DebugClean, t);

if(isnan(DCC_Ind))  %Couldn't find dominant fast frequency, so no DCC
  DCC_Freq = Inf;
else
  [DCC_Freq, f, Pxx, MaxP] = RefineDCC_Freq(t, DV, DCC_Ind, DCC_Corr, ...
					    DebugClean);
end

varargout = {};
if(isfinite(DCC_Freq))
  %Interpolate T, V, and if necessary, I
  T = t(1):(1/DCC_Freq):t(end);
  V = interp1(t, v, T);
  for n = 1:length(varargin)
    varargout{n} = interp1(t, varargin{n}, T);
  end
else
  T = t;
  V = v;
  for n = 1:length(varargin)
    varargout{n} = varargin{n}
  end  
end

if(DeltaArgs == 1)
  DCC_Info.DCC_Freq = DCC_Freq;
  DCC_Info.DCC_Power = MaxP;
  DCC_Info.f = f;
  DCC_Info.Power = Pxx;
  varargout{nargout - 2} = DCC_Info;
end

if(DebugClean)
  h = NamedFigure('Voltage Trace');
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(t, v, 'r.', 'MarkerSize', 6)
  hold on
  plot(T, V, 'b.', 'MarkerSize', 6)
  hold off
  keyboard
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AutoCorr, DV] = GetAutoCorr(t, v, NumAutoCorr)
DV = abs(diff(v));
%DV = abs(PolyDeriv(v, t(2) - t(1)));
[Max1, Ind1] = max(DV);
[MaxDV, TraceInd] = max(Max1);

DV = zscore(DV(:,TraceInd));
MaxLen = 2^20 - 1;
if(length(DV) > MaxLen)
  DV = DV((end-MaxLen+1):end);
end
%The places with big changes are the most useful:
DiscardSmall = false;
DiscardSmall2 = false;
if(DiscardSmall)
  Ind = find(DV > 3);
  if(length(Ind > 1000))
    DV2 = DV;
    DV(:) = 0;
    DV(Ind) = DV2(Ind);
  else
    Ind = find(DV > 2);
    if(length(Ind > 1000))
      DV2 = DV;
      DV(:) = 0;
      DV(Ind) = DV2(Ind);
    else
      Ind = find(DV > 1);
      if(length(Ind > 1000))
	DV2 = DV;
	DV(:) = 0;
	DV(Ind) = DV2(Ind);
      end
    end
  end
end
if(DiscardSmall2)
  Ind = find(DV > 0);
  DV2 = DV;
  DV(:) = 0;
  DV(Ind) = DV2(Ind);
end

AutoCorr = xcorr(DV, NumAutoCorr - 1, 'unbiased');
AutoCorr = AutoCorr(NumAutoCorr:end);
AutoCorr = AutoCorr / AutoCorr(1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DCC_Ind, DCC_Corr] = GetAutoCorrMax(AutoCorr, DebugClean, t)
FSample = 1.0 / (t(2) - t(1));
[Pxx_Corr, f_Corr] = pmtm(AutoCorr, 7/2, [], FSample);
Ind1 = find(f_Corr >= 0.25, 1);
Ind2 = find(f_Corr > 10, 1) - 1;
if(length(Ind2) == 0)
  Ind2 = length(f_Corr);
end
%Pxx_Corr = Pxx_Corr(Ind1:Ind2);
%f_Corr = f_Corr(Ind1:Ind2);

[MaxPow, MaxInd] = max(Pxx_Corr(Ind1:Ind2));
MaxInd = MaxInd + Ind1 - 1;
DCC_Ind = 1 + round(FSample / f_Corr(MaxInd));
if(DCC_Ind == 1)
  DCC_Ind = 2;
elseif(DCC_Ind == length(AutoCorr))
  DCC_Ind = DCC_Ind - 1;
end
if(AutoCorr(DCC_Ind + 1) > AutoCorr(DCC_Ind))
  DCC_Ind = DCC_Ind + 1;
elseif(AutoCorr(DCC_Ind - 1) > AutoCorr(DCC_Ind))
  DCC_Ind = DCC_Ind - 1;
end
DCC_Corr = AutoCorr(DCC_Ind);

%temp for debugging:
Interval = t(DCC_Ind) - t(1);
DCC_Freq = 1.0 / Interval;
%if((DCC_Freq < 0.5 | DCC_Freq > 2.0) & ~(DCC_Corr > 0.3))
%  DebugClean = true;
%end

if(DebugClean)
  h = NamedFigure('CleanDCC Autocorrelation');
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(t(1:length(AutoCorr))-t(1), AutoCorr);
  hold on
  plot(t(DCC_Ind)-t(1), DCC_Corr, 'ro', 'MarkerFaceColor', 'r')
  hold off
  xlabel('Time (ms)', 'FontSize', 18);
  ylabel('Autocorrelation', 'FontSize', 18);
  title('CleanDCC Autocorrelation', 'FontSize', 18);
  
  h = NamedFigure('CleanDCC Autocorrelation Power');
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(f_Corr, Pxx_Corr)
  xlabel('Frequency (kHz)', 'FontSize', 18);
  ylabel('Power', 'FontSize', 18);
  title('CleanDCC Autocorrelation Power', 'FontSize', 18);
end
%Signals to other routines that there is a problem:
if((DCC_Freq < 0.5 | DCC_Freq > 2.0) & ~(DCC_Corr > 0.3))
  DCC_Corr = NaN;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DCC_Freq, f, Pxx, MaxP] = RefineDCC_Freq(t, DV, DCC_Ind, ...
						  DCC_Corr, DebugClean)
if(DCC_Ind <= 2)  %Don't trust these results
  Questionable = true;
  DCC_Freq = 1.0;  %Good a guess as any...
  DCC_Range = 2.0;  %Set a wide range
else
  Interval = t(DCC_Ind) - t(1);
  DCC_Freq = 1.0 / Interval;
  if((DCC_Freq < 0.5 | DCC_Freq > 2.0) & ~(DCC_Corr > 0.3))
    DCC_Freq = 1.0;
    DCC_Range = 2.0;
    Questionable = true;
  else
    Questionable = false;
    DCC_Range = 1.25;
  end
end

if(DebugClean)
  disp(sprintf('Prelim:  Interval = %g, DCC_Freq = %g', Interval, DCC_Freq))
end

FSample = 1.0 / (t(2) - t(1));
try
  [Pxx, f] = pmtm(DV, 7/2, [], FSample);
catch
  DV = DV(1:500000);
  [Pxx, f] = pmtm(DV, 7/2, [], FSample);
end

Ind = find(f > DCC_Freq / DCC_Range & f < DCC_Freq * DCC_Range);

[MaxP, Ind2] = max(Pxx(Ind));
DCC_Freq = f(Ind(Ind2));

if(length(DCC_Freq) == 0)
  disp('Warning:  QUESTIONABLE results in CleanDCC.m')
  DCC_Freq = Inf;
end

if(DebugClean)
  disp(sprintf('DCC Freq = %g kHz', DCC_Freq))
  h = NamedFigure('CleanDCC Power');
  set(h, 'WindowStyle', 'docked');
  plot(f, Pxx);
  hold on
  plot(DCC_Freq, MaxP, 'ro')
  hold off
  xlabel('Frequency (kHz)', 'FontSize', 18)
  ylabel('Power', 'FontSize', 18)
  title('CleanDCC Power', 'FontSize', 18)
end
return