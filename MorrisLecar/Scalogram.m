function varargout = Scalogram(Signal, varargin)
% [FrequencyList, Amplitude, Phase] = Scalogram(Signal, t, PlotSubject)
%
% Calculates the scalogram of Signal and plots it.  The signal is
% required to be sampled at regular time intervals.
%  THEORY:
% Time-Frequency analysis based on a Gabor Transform.
% The basic algorithm is to filter the signal in narrow frequency
% bands.  The width of each band is chosen to be narrow so that the
% time-space width is large-ish (> 1 cycle).  This means that at any
% time, the filtered signal will have a large amplitude if it is
% present, regardless of the original signal's phase.
%  INPUT PARAMETERS:
%   -Signal is 1D vector to be analyzed
%    OPTIONAL:
%    -t is 1D vector of signal times (in seconds)
%    -PlotSubject requests a plot.  Can be boolean true/false
%     (defaults to false), or alternatively it can be a string
%     naming the subject (e.g. 'Waveform 23')
%  OUTPUT PARAMETERS:  (All OPTIONAL)
%   -FrequencyList is list of analyzed frequencies
%   -Amplitude is NumFrequency by NumTime matrix of signal amplitudes
%   -Phase is NumFrequency by NumTime matrix of signal phases


%   Tedious stuff.  Make sure the Signal is properly oriented, and sized.
if(nargin < 1)
  error('Too few arguments.  Run "help Scalogram".')
end
if(size(Signal, 1) > 1)
  if(size(Signal, 2) > 1)
    error('Signal must be 1D vector')
  else
    Signal = Signal';
  end
elseif(size(Signal, 2) <= 1)
  error('Signal must be a 1D vector with more than one element')
end
NumT = length(Signal);

%   More Tedious stuff.  Get all the inputs, check them, etc.
if(nargin == 1)
  t = 0:(NumT-1);
  PlotSubject = '';
  ShowPlot = false;
elseif(nargin == 2)
  if(ischar(varargin{1}))
    PlotSubject = varargin{1};
    ShowPlot = true;
    t = 0:(NumT-1);
  elseif(islogical(varargin{1}))
    PlotSubject = '';
    ShowPlot = varargin{1};
    t = 0:(NumT-1);
  else
    t = varargin{1};
    PlotSubject = '';
    ShowPlot = false;
  end
elseif(nargin == 3)
  t = varargin{1};
  if(ischar(varargin{2}))
    PlotSubject = varargin{2};
    ShowPlot = true;
  else
    PlotSubject = '';
    ShowPlot = varargin{2};
  end
elseif(nargin > 3)
  error('Too many arguments.  Run "help Scalogram"')
end

%  ...Okay.  Now do the Fast Fourier Transform, get frequencies.
SignalFFT = fft(Signal);
DeltaT = t(2) - t(1);
MaxInterval = DeltaT * (NumT - 1);
f = (0:(NumT-1)) / NumT;

% Set some parameters.  These are pretty optimized, so don't mess
% with them.
PlotAmplitude = true & ShowPlot;
PlotPhase = false & ShowPlot;
CalcPhase = (nargout == 3) | PlotPhase;
NumCycles = 2.0; %Time-domain width of filtered frequency, in terms
                 %of number of cycles (MUST be greater than 1, this
                 %is probably best).
MaxFreq = 1.0 / (3 * DeltaT);
MinFreq = 2 * NumCycles / MaxInterval;
NumSigma = 0.5;  %Amount of overlap in neighboring frequency bands

%First, the ideal FreqMult (the spacing between sampled frequencies):
FreqMult = 1.0 + NumSigma / (pi * NumCycles);
NumFreq = ceil(log(MaxFreq / MinFreq) / log(FreqMult)) + 1;
%Next adjust FreqMult to meet demanded MaxFreq and MinFreq:
FreqMult = (MaxFreq/MinFreq).^(1.0/(NumFreq-1));
%Finally, form the list of sampled frequencies
FreqList = (MaxFreq * DeltaT) * FreqMult.^-(0:(NumFreq-1))';

%These help to form the Filter
Filter_Offset = sqrt(.5) * NumCycles * pi;
if(mod(NumT, 2) == 1)
  FLen = (NumT + 1)/2;
  FInd = 1:FLen;
else
  FLen = NumT/2 + 1;
  FInd = 1:FLen;
end
  
f_Scale = f(FInd) * Filter_Offset;

%Pre-allocate for faster speed
Amplitude = zeros(NumFreq, NumT);
if(CalcPhase)
  Phase = zeros(NumFreq, NumT);
end
Filter = zeros(1, NumT);
FilteredSignal = zeros(1, NumT);

%Main loop.  Loop through frequencies
for Row = 1:NumFreq
  Freq = FreqList(Row);
  
  % Construct filter
  Filter_Norm = NumCycles / Freq;  %This is the correct Norm
  %Filter_Norm = 1.0;  %This emphasizes large frequencies (and is wrong)
  Filter(FInd) = Filter_Norm * exp(-(f_Scale / Freq - Filter_Offset).^2);
  
  % Multiply filter and FFT of signal, then take inverse FFT.
  FilteredSig = ifft(SignalFFT .* Filter);
  
  Amplitude(Row,:) = abs(FilteredSig);    % Record the amplitude of the result
  if(CalcPhase)
    Phase(Row,:) = angle(FilteredSig);      % .. and the phase.
  end
end

% Algorithm's all done, now for more tedious stuff!

FreqList = FreqList / DeltaT;
%Pass back the requested results
switch(nargout)
 case 0, varargout = {};
 case 1, varargout = {FreqList};
 case 2, varargout = {FreqList, Amplitude};
 otherwise, varargout = {FreqList, Amplitude, Phase};
end

%Set up tick marks and labels for plots
NYTicks = 8;
dTick = (NumFreq-1)/(NYTicks-1);
YTickLocations = 1:dTick:NumFreq;
YTickValues = interp1(1:NumFreq, FreqList, YTickLocations);

NXTicks = 8;
dTick = (NumT-1) / (NXTicks-1);
XTickLocations = 1:dTick:NumT;
XTickValues = interp1(1:NumT, t, XTickLocations);

%Plotting Amplitude
if(PlotAmplitude)
  NumColors = 256;
  Amplitude = log(Amplitude.*Amplitude);
  MaxVal = max(max(Amplitude));
  MinVal = min(min(Amplitude));
  if(length(PlotSubject) > 0)
    TitleStr = ['Scalogram Amplitude for ', PlotSubject];
  else
    TitleStr = 'Scalogram Amplitude';
  end
  h = NamedFigure(TitleStr);
  set(h, 'WindowStyle', 'docked');
  
  image((NumColors - 1) * (Amplitude - MinVal) / (MaxVal - MinVal));
  RainbowColorMap
  h = gca;
  set(h, 'YTick', YTickLocations);
  set(h, 'YTickLabel', YTickValues);
  set(h, 'XTick', XTickLocations);
  set(h, 'XTickLabel', XTickValues);
  ylabel('Frequency (Hz)', 'FontSize', 18)
  xlabel('Time (s)', 'FontSize', 18)
  title(RealUnderscores(TitleStr), 'FontSize', 18)
end
%Plotting phase
if(PlotPhase)
  if(length(PlotSubject) > 0)
    TitleStr = ['Scalogram Phase for ', PlotSubject];
  else
    TitleStr = 'Scalogram Phase';
  end
  h = NamedFigure(TitleStr);
  set(h, 'WindowStyle', 'docked');
  % Phase only scalogram: phase encoded by hue, saturation uniform.
  hsv(:,:,1) = (Phase+pi)/(2*pi);               % hue varies with phase.
  hsv(:,:,2) = ones(size(Amplitude));           % saturation fixed at 1
  hsv(:,:,3) = ones(size(Amplitude));           % intensity is fixed at 1.
  image(hsv2rgb(hsv));
  h = gca;
  set(h, 'YTick', YTickLocations);
  set(h, 'YTickLabel', YTickValues);
  set(h, 'XTick', XTickLocations);
  set(h, 'XTickLabel', XTickValues);
  ylabel('Frequency (Hz)', 'FontSize', 18)
  xlabel('Time (s)', 'FontSize', 18)
  title(TitleStr, 'FontSize', 18)
end

return
