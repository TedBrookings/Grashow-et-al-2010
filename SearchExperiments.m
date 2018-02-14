function SearchExperiments(varargin)
% SearchExperiments(Flags)
% Looks though analyzed experiments and finds ones that meet a
% criterion, then reports their details and plots the waveforms
%  INPUTS:
%   Flags should be a series of comma separated name, value pairs
%   e.g. SearchExperiments('CellType', 'GM', 'Category', 3, ...
%                          'DutyCycle', [0.25, 0.5])
%        will find all GM cells that are half-center (category 3)
%        with a duty cycle between 0.25 and 0.5
%   available flags:
%    'CellType'  (should be a name, e.g. 'LP')
%    'Category' (0, 1, 2, 3)
%   the rest are ranges (e.g. [0.9, 1.0])
%    'BurstFreq', 'AutoCorr', 'SpikesPerBurst', ...
%    'BurstSpikeFreq', 'DutyCycle', 'SlowWaveAmp'

doPlotWaveforms = false;
close all;

ExperimentsExists = evalin('base', 'exist(''ML_Experiments'')');
if(ExperimentsExists)
  Experiments = evalin('base', 'ML_Experiments');
else
  Experiments = LoadAllExperiments('MorrisLecarFolders.txt');
  assignin('base', 'ML_Experiments', Experiments);
end

try
  MatchList = GetSpecifiedAnalysis(Experiments, varargin{:});
catch errRecord
  help SearchExperiments
  rethrow(errRecord)
end

Intrinsics = LoadIPxls;
ValidID = unique({Intrinsics.ID});

an = MatchList(1);
if(length(an.matchData) == 0)
  matchNames = {};
else
  matchNames = fieldnames(an.matchData);
end

if ispc
  slash = '\';
else
  slash = '/';
end

fprintf('netID\tgsyn\tgh')
for n = 1:length(matchNames)
    fprintf('\t%s', matchNames{n})
end
fprintf('\n')

keepList = [];
for ind = 1:length(MatchList);
  an = MatchList(ind);
  cellID = an.ID;
  if(sum(strcmp(ValidID, cellID)) == 0)
    continue
  else
    keepList = [keepList, ind];
  end
  
  fileName = an.FileName;
  ind1 = strfind(fileName, slash);
  ind1 = ind1(end);
  if ispc
    ind_ = strfind(fileName, '_');
    if(ind_(1) > ind1)
      subdir = [fileName((ind1+1):(ind_(2)-1)), slash];
      fileName = [fileName(1:ind1), subdir, fileName((ind1+1):end)];
      ind1 = strfind(fileName, slash);
      ind1 = ind1(end);
    end
  end
  ind2 = strfind(fileName, '.');
  netID = fileName((ind1+1):(ind2-1));
  
  fprintf('%s\t%d\t%d', netID, an.g_syn, an.g_h)
  for n = 1:length(matchNames)
    fprintf('\t%g', an.matchData.(matchNames{n}))
  end
  fprintf('\n')

  if doPlotWaveforms
    displayWaveform(fileName, netID)
  end
end
MatchList = MatchList(keepList);

CellIDs = unique({MatchList.ID});
fprintf(' Number of unique cells:  %g\n', length(CellIDs))
for n = 1:length(CellIDs)
  fprintf('    %s\n', CellIDs{n})
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayWaveform(fileName, netID)
[t, vReal, vModel] = getData(fileName);

yRange = [min(min(vReal), min(vModel)), max(max(vReal), max(vModel))];

h = NamedFigure(netID);
set(h, 'WindowStyle', 'docked')
clf
h_axes = axes('Position', [0.0 0.0 1.0 1.0], 'Visible', 'off');
xPad = 0.075;
yPad = 0.075;
Width = 1.0 - 2.0 * xPad;
Height = 1.0 - 2.0 * yPad;

x = xPad;
y = yPad + 0.5 * Height;
axes('Position', [x, y, Width, 0.5*Height])  
plot(t, vReal, 'b-', 'LineWidth', 2);
set(gca, 'xTickLabel', {})
ylabel('Real Cell Potential (mV)', 'FontSize', 18')
title(RealUnderscores(netID), 'FontSize', 18')
ylim([-75 -5])
xlim([0 40])

x = xPad;
y = yPad;
axes('Position', [x, y, Width, 0.5*Height])  
plot(t, vModel, 'r-', 'LineWidth', 2);
ylabel('Model Cell Potential (mV)', 'FontSize', 18')
xlabel('Time (s)', 'FontSize', 18')
ylim([-75 -5])
xlim([0 40])
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, vReal, vModel] = getData(fileName)
AbfS = LoadAbf(fileName);
tBig = AbfS.Time; % * 1000;  %convert to ms
t = tBig(1):0.001:tBig(end);

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
  vReal = AbfS.Data.(FNames{Voltage(2)});
  vModel = AbfS.Data.(FNames{Voltage(1)});
elseif(StringCheck(FNames{Voltage(2)}, 'model'))
    vReal = AbfS.Data.(FNames{Voltage(1)});
    vModel = AbfS.Data.(FNames{Voltage(2)});
else
  disp(FNames)
  error(['No model traces found in ', FileName])
end

vReal = interp1(tBig, vReal, t);
vModel = interp1(tBig, vModel, t);
return