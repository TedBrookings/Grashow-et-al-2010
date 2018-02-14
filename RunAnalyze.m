function Analysis = RunAnalyze(FileName, PlotSubject)
% Analysis = RunAnalyze(FileName, PlotSubject)
% (1) loads the .smr file 'FileName',
% (2) examines two voltage traces (presumably from a dynamic clamp
%      experiment) and calculates various properties of each trace
%      and their mutual interactions, and
% (3) categorizes the results
%
%  INPUT PARAMETERS:
%    -FileName is the name (with complete path) of the .smr file
%   OPTIONAL:
%     -PlotSubject should be set to true[false] to produce[suppress]
%       plots of waveforms/analysis.  Alternatively, it can be set
%       to a string to aid it titling plots (e.g. 'Exp #71')
%       PlotVar defaults to false
%  OUTPUT PARAMETERS:
%    -Analysis.Cat:  integer value that codes the category type.
%       type 'help CategorizeDynamic' for info.
%    -Analysis.CatString: string that describes the category
%    -Analysis.Cell0:  structure with information about cell 0
%       returned from AnalyzeWaveform.m.  It should be
%       straightforward, but type 'help AnalyzeWaveform3' for info.
%    -Analysis.Cell1:  structure with information about cell 1
%    -Analysis.SlowWave:  structure with information about
%     slow-wave behaviour of the combined system
%       -SlowWave.Freq is the frequency of the dominant slow-wave
%        component (in Hz)
%       -SlowWave.Sigma is a (very crude) measure of the importance
%        of the slow-wave frequency in the power spectrum
%    -Analysis.HalfCenter: structure with information about
%     half-center oscillations
%       -HalfCenter.Freq is the half-center frequency (in Hz)
%       -HalfCenter.Phase is the phase that cell 1 lags behind cell0
%       -HalfCenter.PhaseSigma is the standard deviation of the phases
%       -HalfCenter.ExclusionFact is the tendancy of bursts to
%        avoid each other.  0 is random, 1 is maximum, and -1 is minimum.
%
% If a feature is not detected, relevant frequencies are set to
% zero, and relevant lists are empty


StabilizationTime = 10; %s
LowThresh = -1;  %Threshold for spike derivatives, mv/ms
HighThresh = 2;  %  Currently these are ignored!
if(nargin < 2)
  PlotSubject = false;
elseif(~ischar(PlotSubject))
  if(PlotSubject)
    if(ispc)
      Slash = '\';
    else
      Slash = '/';
    end
    Ind = strfind(FileName, Slash);
    Ind = Ind(end) + 1;
    PlotSubject = FileName(Ind:end);
  end
end

[t, v0, v1] = GetData(FileName);

if(DoPlot(PlotSubject))
  if(ischar(PlotSubject) & length(PlotSubject) > 0)
    TitleStr = [PlotSubject, ' Waveforms'];
  else
    TitleStr = 'Waveforms';
  end
  h = NamedFigure(TitleStr);
  set(h, 'WindowStyle', 'docked');
  hold off;
  plot(t / 1000, v0, 'b-');
  hold on;
  plot(t / 1000, v1, 'r-');
  xlabel('Time (s)', 'FontSize', 18);
  ylabel('Voltage (mV)', 'FontSize', 18);
  title(RealUnderscores(TitleStr), 'FontSize', 18);
  hold off;
end

Analysis = AnalyzeDynamic(t, v0, v1, StabilizationTime, ...
			  LowThresh, HighThresh, PlotSubject);
[Cat, CatString] = CategorizeDynamic(Analysis);
disp(sprintf('System is %s', CatString))

Analysis.Cat = Cat;
Analysis.CatString = CatString;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, v0, v1] = GetData(FileName)
if(length(strfind(lower(FileName), '.abf')) > 0)
  [t, v0, v1] = GetCleanAbfData(FileName);
elseif(length(strfind(lower(FileName), '.smr')) > 0)
  [t, v0, v1] = GetCleanSmrData(FileName);
elseif(length(strfind(lower(FileName), '.dat')) > 0)
  fid = fopen(FileName);
  if(fid < 0)
    error(['Couldn''t open ', FileName])
  end
  try
    NumT = fscanf(fid, '%d', 1);
    NumCol = 3;
    Mat = fscanf(fid, '%g', [NumCol, NumT]);
    fclose(fid);
    t = Mat(1,:);
    v0 = Mat(2,:);
    v1 = Mat(3,:);
  catch
    ErrStruct = lasterr;
    disp(ErrStruct.message)
    error(['Error reading ', FileName])
  end
elseif(length(strfind(lower(FileName), '.mat')) > 0)
  load(FileName, 't', 'v0', 'v1');
else
  ErrStr = sprintf(['Error opening ', FileName, ... 
		    '\nFiles must be of type .abf .smr .dat or .mat']);
  error(ErrStr);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, v0, v1] = GetCleanAbfData(FileName)
AbfS = LoadAbf(FileName);
t = AbfS.Time * 1000;  %convert to ms

FNames = fieldnames(AbfS.Units);
Current = [];
Voltage = [];
for n = 1:length(FNames)
  Unit = AbfS.Units.(FNames{n});
  if(strcmp(Unit, 'mV'))
    Voltage = [Voltage, n];
  elseif(strcmp(Unit, 'nA'))
    Current = [Current, n];
  end
end
if(length(Voltage == 2) ~= 2)
  for n = 1:length(FNames)
    disp(FNames)
  end
  error(['Incorrect number of voltage traces in ', FileName]);
end
  
v0 = AbfS.Data.(FNames{Voltage(1)});
v1 = AbfS.Data.(FNames{Voltage(2)});

[t, v0, v1] = CleanAndSmooth(t, v0, v1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, v0, v1] = GetCleanSmrData(FileName)
[t, v_arr, ChanNames] = load_smr(FileName);

%Get index of voltage channels
Ind = [];
NameInd = strfind(ChanNames, '1V');
for n = 1:length(NameInd)
  if(length(NameInd{n}) > 0)
    Ind = [Ind, n];
  end
end
if(length(Ind) ~= 2)
  ChanNames
  error(['Incorrect number of voltage traces in ', FileName]);
end

if(size(v_arr, 1) ~= length(t))
  v_arr = v_arr';
end

v0 = v_arr(:,Ind(1));
v1 = v_arr(:,Ind(2));
%t = t * 1000;
while(t(2) - t(1) < .01)  %this is a kludge!
  t = t * 1000;
end

%This step necessary with one-electrode set-up:
[t, v0, v1] = CleanAndSmooth(t, v0, v1);
return
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, v0, v1] = CleanAndSmooth(t, v0, v1)
%First remove DCC noise by interpolating to the DCC frequency
[t, v0, v1, DCC_Info] = CleanDCC(t, v0, v1);
if(~isfinite(DCC_Info.DCC_Freq))
  disp('Warning:  weird DCC signal!')
  DCC_Info
end
disp(sprintf('DCC freq = %g kHz', DCC_Info.DCC_Freq))

%Next low-pass filter the waveform
SmoothTime = .5;  %ms
if(t(2) - t(1) < SmoothTime)
  n = 2;
  while(t(n) - t(1) < SmoothTime)
    n = n + 2;
  end
  [B,A] = butter(2, 2 / n,'low');
  v0 = filtfilt(B, A, v0);
  v1 = filtfilt(B, A, v1);
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