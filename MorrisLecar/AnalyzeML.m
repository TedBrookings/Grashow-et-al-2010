function Analysis = AnalyzeML(varargin)
% Analysis = AnalyzeML(FileName, PlotSubject)
%       ...OR...
% Analysis = AnalyzeML(t, VReal, VModel, PlotSubject)
%
% NOTE: Second form is deprecated.
% Examines two voltage traces, which are assumed to be from a
% dynamic clamp experiment where the first cell is a real neuron,
% and the second is a Morris-Lecar model cell.
% AnalyzeML calculates various properties of each trace and
% properties of their mutual interactions.
%
%  INPUT PARAMETERS:
%   -t is time in ms
%   -VReal is voltage of real cell (in mV)
%   -VModel is voltage of Morris-Lecar model cell
%    OPTIONAL:
%     -PlotSubject should be set to true[false] to produce[suppress]
%       plots of waveforms/analysis.  Alternatively, it can be set
%       to a string to aid it titling plots (e.g. 'Exp #71')
%       PlotVar defaults to false
%  OUTPUT PARAMETERS:
%    -Analysis.CellReal:  structure with information about real cell
%       returned from AnalyzeWaveform3.m.  It should be
%       straightforward, but type 'help AnalyzeWaveform3' for info.
%    -Analysis.SlowWave:  structure with information about
%     slow-wave behaviour of the combined system
%       -SlowWave.Freq: frequency of the dominant slow-wave
%        component (in Hz)
%       -SlowWave.Sigma: (very crude) measure of the importance
%        of the slow-wave frequency in the power spectrum
%       -SlowWave.Corr: autocorrelation at slow-wave period
%       -SlowWave.Scalogram: structure with spectrum information
%         -Scalogram.FrequencyList:  list of analyzed frequencies
%         -Scalogram.AvgPower: length NumFreq list of average powers
%         -Scalogram.SlowPower: length NumFreq list of powers
%                sampled only between spikes
%         NOTE:  The next two are not returned if a plot is not
%         requested, because they take up a LOT of space:
%         -Scalogram.Amplitude:  NumFreq x NumTime matrix of
%             scalogram amplitudes
%         -Scalogram.Phase:  NumFreq x NumTime matrix of phases
%
% If a feature is not detected, relevant frequencies are set to
% zero, and relevant lists are empty

[t, VReal, VModel, SynapseData, PlotVar, PlotSubject] = GetArgs(varargin);

%close all;
PlotBrief = true;  %add text(x,y, 'string');
if(PlotVar)
  if(length(PlotSubject) > 0)
    TitleStr = [PlotSubject, ' Waveforms'];
  else
    TitleStr = 'Waveforms';
  end
  MinV = min([min(VReal), min(VModel)]);
  MaxV = max([max(VReal), max(VModel)]);
  h = NamedFigure(TitleStr);
  hold off
  set(h, 'WindowStyle', 'docked');
  subplot('Position', [0.05, 0.51, 0.945, 0.42])
  plot(t/1000, VReal, 'b-')
  ylim([MinV, MaxV])
  set(gca, 'xTickLabel', {})
  title(RealUnderscores(TitleStr), 'FontName', 'Ariel', 'FontSize', 18)
  %xlabel('Time (s)', 'FontSize', 18)
  ylabel('Voltage (mV)', 'FontName', 'Ariel', 'FontSize', 18)
  hold on
  subplot('Position', [0.05, 0.09, 0.945, 0.42])
  plot(t/1000, VModel, 'r-')
  ylim([MinV, MaxV])
  xlabel('Time (s)', 'FontName', 'Ariel', 'FontSize', 18)
  ylabel('Voltage (mV)', 'FontName', 'Ariel', 'FontSize', 18)
  hold off
  
  if(length(PlotSubject) > 0)
  TitleStr = [PlotSubject, ' Synaptic Current'];
  else
    TitleStr = 'Synaptic Current';
  end
  h = NamedFigure(TitleStr);
  hold off
  set(h, 'WindowStyle', 'docked');  
  plot(t/1000, SynapseData.I_FromModelToBio, 'b-')
  hold on
  plot(t/1000, SynapseData.I_FromBioToModel, 'r-')
  hold off
  xlabel('Time (s)', 'FontName', 'Ariel', 'FontSize', 18)
  ylabel('Current (nA)', 'FontName', 'Ariel', 'FontSize', 18)  
  title(RealUnderscores(TitleStr), 'FontName', 'Ariel', 'FontSize', 18)
end

StabilizationTime = 10; %s

nStable = find(t > t(1) + StabilizationTime * 1000, 1);
t = t(nStable:end);
VReal = VReal(nStable:end);
VModel = VModel(nStable:end);
SynapseData.I_FromModelToBio = SynapseData.I_FromModelToBio(nStable:end); 
SynapseData.I_FromBioToModel = SynapseData.I_FromBioToModel(nStable:end);

CellReal = AnalyzeWaveform(t, VReal, PlotSubject);
if(PlotVar & ischar(PlotSubject) & length(PlotSubject) > 0)
  ModelSubject = [PlotSubject, ' Model'];
else
  ModelSubject = PlotSubject;
end
ModelSlow = AnalyzeSlowWave(t, VModel, [], ModelSubject);
ModelSlow.Phases = [];  %Reduce storage demand

%Check to make sure that the "bursts" are halfcenter
CellReal = CheckHalfCenter(CellReal, t, VReal, VModel, SynapseData);

%Now this should be calculated by AnalyzeWaveform.m
disp(sprintf('%s = %g Hz, %s = %g Hz, %s = %g Hz', ...
	     'SpikeFreq', CellReal.Spike.Freq, ...
	     'BurstFreq', CellReal.Burst.Freq, ...
	     'BurstSpikeFreq', CellReal.Burst.SpikeFrequencies.Mean))
disp(sprintf('SlowFreq = %g Hz, Autocorrelation = %g', ...
	     CellReal.SlowWave.Freq, CellReal.SlowWave.Corr))
disp(sprintf('Model SlowFreq = %g Hz, Autocorrelation = %g', ...
	     ModelSlow.Freq, ModelSlow.Corr))

Analysis.CellReal = CellReal;
Analysis.ModelSlow = ModelSlow;
[Cat, CatString] = CategorizeML(Analysis);
disp(sprintf('System is %s', CatString))

Analysis.Cat = Cat;
Analysis.CatString = CatString;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, VReal, VModel, SynapseData, PlotVar, PlotSubject] = GetArgs(Args)

NumArgs = length(Args);
if(NumArgs < 1)
  error('Too few arguments.  Run "help AnalyzeML".')
end
if(ispc)
  Slash = '\';
else
  Slash = '/';
end
PlotSubject = '';
if(ischar(Args{1}))  %passed a filename
  FileName = Args{1};
  [t, VReal, VModel, SynapseData] = GetData(FileName);
  if(NumArgs == 1)
    PlotVar = false;
  elseif(NumArgs == 2)
    PlotVar = Args{2};
    Ind = strfind(FileName, Slash);
    Ind = Ind(end) + 1;
    PlotSubject = FileName(Ind:end);
  else
    error('Passed too many arguments.  Run "help AnalyzeML".')
  end
else   %passed data
  if(NumArgs < 3)
    error('Invalid arguments.  Run "help AnalyzeML".')
  end
  t = Args{1};
  VReal = Args{2};
  VModel = Args{3};
  if(NumArgs == 3)
    PlotVar = false;
  elseif(NumArgs == 4)
    PlotVar = Args{4};
  else
    error('Passed too many arguments.  Run "help AnalyzeML".')
  end
  fprintf(2, 'Warning, passing data to AnalyzeML is deprecated.\n')
  g_syn = 10;
  SynapseData = SimSynapses(t, VReal, VModel, g_syn);
end
if(PlotVar == false)
  PlotSubject = false;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, VReal, VModel, SynapseData] = GetData(FileName)
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
  error(['Incorrect number of voltage traces in ', FileName])
end

if(StringCheck(FNames{Voltage(1)}, 'model'))
  if(StringCheck(FNames{Voltage(2)}, 'model'))
    disp(FNames)
    error(['Two model traces found in ', FileName])
  end
  VReal = AbfS.Data.(FNames{Voltage(2)});
  VModel = AbfS.Data.(FNames{Voltage(1)});
elseif(StringCheck(FNames{Voltage(2)}, 'model'))
    VReal = AbfS.Data.(FNames{Voltage(1)});
    VModel = AbfS.Data.(FNames{Voltage(2)});
else
  disp(FNames)
  error(['No model traces found in ', FileName])
end

g_syn = Get_g_syn(FileName);
SynapseData = SimSynapses(t, VReal, VModel, g_syn);
[t, VReal, VModel, SynapseData] = CleanAndSmooth(t, VReal, VModel, SynapseData);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g_syn = Get_g_syn(FileName)
if(ispc)
  Slash = '\';
else
  Slash = '/';
end
IndSlash = strfind(FileName, Slash);
DirName = FileName(1:IndSlash(end));

FileList = dir(DirName);
FileList = GetTextFiles(FileList);

if(length(FileList) ~= 1)
  error('Directory must have exactly one .txt file.');
end

ExpFileName = [DirName, FileList.name];


IndUnder = strfind(FileName, '_');
Cat = FileName(IndUnder(3)+1:IndUnder(4)-1);
IndDot = strfind(FileName, '.');
FileNum = str2num(FileName(IndUnder(4)+1:IndDot(1)-1));

fid = fopen(ExpFileName);
if(fid < 0)
  error(['Error opening ', ExpFileName]);
end

dummy = fgets(fid);  %contains the header description

g_syn = NaN;
while(true)
  Temp = fscanf(fid, '%d', 3);
  if(length(Temp) < 3)
    break;
  end

  ThisFileNum = Temp(3);
  ThisCat = fscanf(fid, '%s', 1);
  
  if(ThisFileNum ~= FileNum || ~strcmp(ThisCat, Cat))
    continue;
  end
  
  g_syn = Temp(1);
end
fclose(fid);
if(isnan(g_syn))
  fprintf(2, 'Warning:  Couldn''t find %s File %g in %s\n', ...
              Cat, FileNum, ExpFileName)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TextFileList = GetTextFiles(FileList)
TextFileList = [];
for n = 1:length(FileList)
  if(strfind(FileList(n).name, '.txt'))
    TextFileList = [TextFileList, FileList(n)];
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, VReal, VModel, SynapseData] = CleanAndSmooth(t, VReal, VModel, ...
                                                          SynapseData)
%First remove DCC noise by interpolating to the DCC frequency
[t, VReal, VModel, ...
  SynapseData.I_FromBioToModel, SynapseData.I_FromModelToBio, ...
  DCC_Info] = CleanDCC(t, VReal, VModel, ...
  SynapseData.I_FromBioToModel, SynapseData.I_FromModelToBio);

if(~isfinite(DCC_Info.DCC_Freq))
  disp('Warning:  weird DCC signal!')
  DCC_Info
end
fprintf('DCC freq = %g kHz\n', DCC_Info.DCC_Freq)

%Next low-pass filter the waveform
SmoothTime = .5;  %ms
if(t(2) - t(1) < SmoothTime)
  n = 2;
  while(t(n) - t(1) < SmoothTime)
    n = n + 2;
  end
  %consider changing butter to besself, probably of higher order
  [B,A] = butter(2, 2 / n,'low');
  VReal = filtfilt(B, A, VReal);
  VModel = filtfilt(B, A, VModel);
end
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
function CellReal = CheckHalfCenter(CellReal, t, VReal, VModel, SynapseData)
if(CellReal.Burst.Freq <= 0)  %No burst, no half center
  CellReal.HalfCenter = false;
  return
end

nOn = CellReal.Burst.StartInds;
nOff = CellReal.Burst.StopInds;

ExclusionFrac = 0.9;
HigherInBurst = 0;
TotalInBurst = 0;
HigherOutBurst = 0;
TotalOutBurst = 0;
%Note, < comparisons are correct, because the current is negative
for m=1:(length(nOn)-1)
  Ind = nOn(m):nOff(m);
  TotalInBurst = TotalInBurst + 1;
  Ind = find(SynapseData.I_FromBioToModel(Ind) < ...
             SynapseData.I_FromModelToBio(Ind));
  HigherInBurst = HigherInBurst + (length(Ind) > 0);

  Ind = nOff(m):nOn(m+1);
  TotalOutBurst = TotalOutBurst + 1;
  Ind = find(SynapseData.I_FromModelToBio(Ind) < ...
             SynapseData.I_FromBioToModel(Ind));  
  HigherOutBurst = HigherOutBurst + (length(Ind) > 0);
end
m = length(nOn);
Ind = nOn(m):nOff(m);
TotalInBurst = TotalInBurst + 1;
Ind = find(SynapseData.I_FromBioToModel(Ind) < ...
           SynapseData.I_FromModelToBio(Ind));
HigherInBurst = HigherInBurst + (length(Ind) > 0);

%FracInBurst = fraction of bursts when synaptic current from biological cell
%              is momentarily greater than synaptic current from model
FracInBurst = HigherInBurst / TotalInBurst;
%FracOutBurst = fraction of bursts when synaptic current from model
%              is momentarily greater than synaptic current from bio cell
FracOutBurst = HigherOutBurst / TotalOutBurst;
fprintf('InBurst: %g, OutBurst: %g\n', FracInBurst, FracOutBurst)

if(FracInBurst < ExclusionFrac || FracOutBurst < ExclusionFrac)
  CellReal.HalfCenter = false;
else
  CellReal.HalfCenter = true;
end

return