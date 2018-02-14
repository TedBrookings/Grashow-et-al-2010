function PlotPooled(varargin)
PlotType = 'ML';
switch(PlotType)
 case 'GMGM', PlotPooledGMGM(varargin{:})
 case 'ML', PlotPooledML(varargin{:})
 otherwise, error(['Invalid plot type ', PlotType])
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPooledGMGM(Organized)

ValueString = 'HalfCenter.Freq';  %Tells which values to map
if(nargin < 1)
  OrganizedExists = evalin('base', 'exist(''Organized'')');
  if(OrganizedExists)
    Organized = evalin('base', 'Organized');
  else
    Organized = OrganizeByParameters([], ValueString);
    assignin('base', 'Organized', Organized);
  end
end

close all

DrawMeanChange(Organized);
DrawScatter(Organized);
DrawBurstFraction(Organized, 'GMGM');
DrawBar(Organized);

Props = LoadIntrinsicProperties;
DrawCorrelateProps(Organized, Props);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPooledML(Experiments)
if(nargin < 1)
  ExperimentsExists = evalin('base', 'exist(''ML_Experiments'')');
  if(ExperimentsExists)
    Experiments = evalin('base', 'ML_Experiments');
  else
    Experiments = LoadAllExperiments('MorrisLecarFolders.txt');
    assignin('base', 'ML_Experiments', Experiments);
  end
end

close all

%///////////// Select from lots of options of what to do ////////////////

%options.netType = 1;   %Cell spiking
options.netType = 3;  %Bursting
%options.netType = 0.85;  %Highly-regular slow wave

options.numSynDivisions = 2;

%Select which kinds of analysis to do:
options.tryIDCellTypes = true;
options.tryCorrelations = false;

%Select which kinds of plots to make:
options.plotIndividual = false;
options.plotImportance = false;
options.plotBurstFraction = false;

%Choose the number of shuffled trials (should be large or zero):
options.numShuffledTrials = 1000;

%Choose which scoreHandle to use:
options.scoreHandle = @ZScore;
%options.scoreHandle = @RankScore;

%Choose which map properties to correlate with IPs/pyloric props
options.corrProps = { ...
    %'Mean f', 'Std f', ...
%    
    %'Frac HalfCenter', ...
    %'Mean g_syn', 'Mean g_h', 'Mean f', ...
    %'Std g_syn', 'Std g_h', 'Std f', ...
    
    %'g_s_y_n Sensitivity', 'g_h Sensitivity', ...
%    
    %DenseLabels{:}, ...
%    
    'Mean AutoCorr', 'Mean SpikesPerBurst', ...
    'Mean BurstSpikeFreq', 'Mean DutyCycle', ...
    'Mean SlowWaveAmp', ...
%    
    %'Std AutoCorr', 'Std SpikesPerBurst', ...
    %'Std BurstSpikeFreq', 'Std DutyCycle', ...
    %'Std SlowWaveAmp' ...
		    };
%Choose which map properties to use to ID cell types:
options.clustProps = ...
    {'Frac HalfCenter', 'Mean f', 'Span g_h', 'Mean SpikesPerBurst'};    
options.netLabels = {'BurstFrequency', 'DutyCycle', 'Phase', ...
		    'SpikesPerBurst'};

%Choose which pyloric properties to use for correlations:

%////////////////////// Load in all the data //////////////////////////
%Load intrinsic properties
intrinsics = LoadIPxls;

%Get lots of properties from the maps:
[mapProps, mapPropLabels] = getAllMapProps(Experiments, options);
%Get Pyloric props for LP cells
LPNetStats = LoadNetStats_xls('LP');

%Get Pyloric props for PD cells
PDNetStats = LoadNetStats_xls('PD');

%Get a list of maps (one map per cell) of networks that meet the
%  netType criterion
mapList = getMapList(Experiments, intrinsics, options);

%///////////////////// Organize data into structures ////////////////////////
%  Among other things, these structures have their data organized
%   so that their rows all correspond to the same individual cells.
IPs = makeIPStruct(intrinsics, options);

corrMapProps = makeMapPropStruct(mapProps, mapPropLabels, intrinsics, ...
				 options, options.corrProps);
sepMapProps = makeMapPropStruct(mapProps, mapPropLabels, intrinsics, ...
				options, options.clustProps);

netProps = makePyloricStruct(LPNetStats, PDNetStats, ...
			     intrinsics, options, ...
			     options.netLabels);

%/////////////////////// Actually do stuff ////////////////////////////

%Orphaned Code, ANOVA's on IPs
%{
  %Do some ANOVAs for each intrinsic property
  cellTypes = {Intrinsics.BaseCond};
  for n = 1:length(XLabels)
    [p,tbl,stats] = anova1(X(:,n), cellTypes);
    %p = p * length(XLabels);
    fprintf('ANOVA %s, p=%g\n', XLabels{n}, p)
    if(p < 0.05)
      [c, m, h, gNames] = multcompare(stats);
      for m=1:size(c,1)
	row = c(m,:);
	Cell1 = gNames{row(1)};
	Cell2 = gNames{row(2)};
	if(row(3) > 0 || row(5) < 0)
	  %Cell1 and Cell2 are signiicantly different for this IP
	  fprintf('\t%s and %s are significantly different\n', ...
		  Cell1, Cell2)
	end
      end
    end
  end

%}


if(options.plotBurstFraction)
  DrawBurstFraction(Experiments, 'ML', intrinsics);
end

if(options.tryIDCellTypes)
  %Try to figure out the cell types from intrinsic/map properties
  SeparateCellTypes(IPs, sepMapProps, mapList, options);
end

if(options.tryCorrelations)
  %Find correlations between intrinsic properties and map properties,
  % then assess significance using shuffled trials, if
  % numShuffledTrials > 0
  runShuffledCorr(options, IPs, corrMapProps);

  %Find correlations between intrinsic properties and network pyloric
  % properties, then assess significance using shuffled trials, if
  % numShuffledTrials > 0
  options.fitPooled = false;
  options.fitSubTypes = true;
  runShuffledCorr(options, IPs, netProps);

  %Find correlations between map properties and network pyloric
  % properties, then assess significance using shuffled trials, if
  % numShuffledTrials > 0
  runShuffledCorr(options, corrMapProps, netProps);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mapProps, mapPropLabels] = getAllMapProps(Experiments, options)
[synSens, hSens] = GetSensitivity(Experiments, options.netType);
SensLabels = {'g_s_y_n Sensitivity', 'g_h Sensitivity'};
[LocalDensities, DenseLabels] = GetLocalDensities(Experiments, ...
						  options.netType, ...
						  options.numSynDivisions);
[MapQuantities, MapQuantitiesLabels] = ...
    GetMapQuantities(Experiments, options.netType);

%Put all the properties together:
mapProps = {MapQuantities{:}, synSens, hSens, LocalDensities{:}};
mapPropLabels = {MapQuantitiesLabels{:}, SensLabels{:}, DenseLabels{:}};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind = getWantedInd(labels, wantedLabels)
ind = [];
for n = 1:length(wantedLabels)
  k = find(strcmp(labels, wantedLabels{n}));
  if(length(k) ~= 1)
    fprintf(2, 'Searching through labels:\n')
    for m = 1:length(labels)
      fprintf(2, '\t%s\n', labels{m})
    end
    error('Error finding %s in labels.', wantedLabels{k})
  end
  ind = [ind, k];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mapList = getMapList(Experiments, intrinsics, options)
numCells = length(intrinsics);
IDList = {intrinsics.ID};
if(options.netType == round(options.netType))
  %tell GetSpecifiedAnalysis to get networks of category options.netType
  flag = 'Category';
  criterion = options.netType;
else
  %tell GetSpecifiedAnalysis to get networks with AutoCorr between
  %  options.netType and Inf
  flag = 'AutoCorr';
  criterion = [options.netType, Inf];
end

mapLen = 7;  %length of map along each dimension
gMin = 10;  %starting value for gsyn and gh
gStep = 15;  %difference in values for neighboring gsyn, gh

mapList = repmat(LocationMap(zeros(mapLen, mapLen)), numCells, 1);
burstList = GetSpecifiedAnalysis(Experiments, flag, criterion);

for n = 1:length(burstList)
  burstStruct = burstList(n);
  mapNum = find(strcmp(burstStruct.ID, IDList));
  if(length(mapNum) == 0)
    continue
  elseif(length(mapNum) > 1)
    error('Multiple cells match ID: %s', burstStruct.ID)
  end
  ind_syn = 1 + round((burstStruct.g_syn - gMin) / gStep);
  ind_h = 1 + round((burstStruct.g_h - gMin) / gStep);
  
  if(ind_syn < 1 || ind_syn > mapLen || ind_h < 1 || ind_h > mapLen)
    fprintf(2, 'Strange index for ID = %s, (gsyn=%g,gh=%g)\n', ...
	    burstStruct.ID, burstStruct.g_syn, burstStruct.g_h)
    continue
  end
  mapList(mapNum).mapMatrix(ind_syn, ind_h) = 1;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IPStruct = makeIPStruct(intrinsics, options)
%Organize intrinsic properties into a matrix,
%  IPStruct.mat which is numCells x numIPs
IP_R = options.scoreHandle(cat(1, intrinsics.R));
IP_V = options.scoreHandle(cat(1, intrinsics.VThresh));
IP_F = options.scoreHandle(cat(1, intrinsics.FISlope));
IP_I = options.scoreHandle(cat(1, intrinsics.SpikeRate1nA));
IP_VRest = options.scoreHandle(cat(1, intrinsics.VRest));
IP_SH = options.scoreHandle(cat(1, intrinsics.SpikeHeight));

IPStruct.label = 'Intrinsic Props';
IPStruct.mat = [IP_R, IP_V, IP_F, IP_I, IP_VRest, IP_SH];
IPStruct.labels = {'R', 'VThresh', 'FI Slope', 'Spike Rate 1nA', ...
		   'VMin', 'Spike Height'};
IPStruct.cellType = {intrinsics.BaseCond};
IPStruct.ID = {intrinsics.ID};
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mapPropStruct = makeMapPropStruct(mapProps, mapPropLabels, ...
					   intrinsics, options, wantedProps)
%This function creates a matrix whos rows represent cells, and
%  correspond to the same rows in the intrinsics matrix

%First select only the properties we care about
wantedInds = getWantedInd(mapPropLabels, wantedProps);
mapProps = mapProps(wantedInds);
mapPropLabels = mapPropLabels(wantedInds);

%Now construct the matrix, leaving NaNs where information
%  doesn't exist
numCells = length(intrinsics);
numMapProps = length(mapProps);
propMat = repmat(NaN, numCells, numMapProps);
for n = 1:numCells
  ID = intrinsics(n).ID;
  for m = 1:numMapProps
    propMatInd = find(strcmp({mapProps{m}.ID}, ID));
    if(length(propMatInd) == 0)
      fprintf(2, 'Warning: No %s data for cell ID: %s\n', ...
	      mapPropLabels{m}, ID)
      %propMat(n,m) will remain NaN
    elseif(length(propMatInd) > 1)
      %this should never happen (?)
      error('Multiple %s data for cell ID: %s', mapPropLabels{m}, ID)
    else
      propMat(n,m) = mapProps{m}(propMatInd).Mean;
    end
  end
end

mapPropStruct.label = 'Map Props';
%mapPropStruct.mat = options.scoreHandle(propMat);
mapPropStruct.mat = propMat;
mapPropStruct.labels = wantedProps;
mapPropStruct.cellType = {intrinsics.BaseCond};
mapPropStruct.ID = {intrinsics.ID};
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pyloricStruct = makePyloricStruct(LPNetStats, PDNetStats, ...
					   intrinsics, options, ...
					   commonLabels)
numCells = length(intrinsics);
numNet = length(commonLabels);
netProps = repmat(NaN, numCells, 2 * numNet);
for n = 1:numCells
  ID = intrinsics(n).ID;
  mLP = find(strcmp(ID, {LPNetStats.ID}));
  mPD = find(strcmp(ID, {PDNetStats.ID}));
  if length(mLP) + length(mPD) == 0
    if ismember(intrinsics(n).BaseCond, {'LP', 'PD'})
      fprintf(2, 'Warning: No pyloric data for ID: %s\n', ID)
    end
  elseif length(mLP) + length(mPD) > 1
    error('Multiple pyloric data for ID: %s\n', ID);
  elseif length(mLP) > 0
    for k = 1:length(commonLabels)
      netProps(n, k) = LPNetStats(mLP).(commonLabels{k});
    end
  else
    for k = 1:length(commonLabels)
      netProps(n, numNet + k) = PDNetStats(mPD).(commonLabels{k});
    end    
  end
end

pyloricStruct.label = 'Pyloric Props';
pyloricStruct.mat = options.scoreHandle(netProps);
netLabels = cell(2 * numNet, 1);
for k = 1:length(commonLabels)
  netLabels{k} = ['LP ', commonLabels{k}];
  netLabels{numNet + k} = ['PD ', commonLabels{k}];
end
pyloricStruct.labels = netLabels;
pyloricStruct.cellType = {intrinsics.BaseCond};
pyloricStruct.ID = {intrinsics.ID};
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = runShuffledCorr(options, xData, yData)

%Do the analysis with the selected map properties
resultStruct = CorrelateProperties(xData, yData, options);

if(options.numShuffledTrials > 0)
  %Make a popup to display progress
  popupLabel = sprintf('%s vs %s shuffled trials', ...
		       xData.label, yData.label);
  PopupProgress(popupLabel, options.numShuffledTrials)
  
  resultStruct.shufflePVals = zeros(resultStruct.numFits, ...
				    options.numShuffledTrials);
  resultStruct.shuffleRSquared = zeros(resultStruct.numFits, ...
				       options.numShuffledTrials);  
  fprintf('\nRunning shuffled trials ...')
  for n = 1:options.numShuffledTrials
    resultStruct.shuffleNum = n;
    resultStruct = CorrelateProperties(xData, yData, options, ...
				       resultStruct);
    PopupProgress(popupLabel)
  end
  f = findobj('Name', popupLabel);
  close(f)
  fprintf(' done.\n')
else
  resultStruct.shufflePVals = [];
  resultStruct.shuffleRSquared = [];
end

CorrectAndDisplay(resultStruct, options)

if nargout == 0
  varargout = {};
else
  varargout = {resultStruct};
end
return