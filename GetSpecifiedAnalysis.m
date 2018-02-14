function MatchList = GetSpecifiedAnalysis(Experiments, varargin)
% MatchList = GetSpecifiedAnalysis(Experiments, Flags)
% Looks though analyzed Experiments and finds ones that meet a
% criterion, then returns them as MatchList
%  INPUTS:
%   Experiments: an array of structures
%   Flags: series of comma separated name, value pairs
%     e.g. SearchExperiments('CellType', 'GM', 'Category', 3, ...
%                          'DutyCycle', [0.25, 0.5])
%        will find all GM cells that are half-center (category 3)
%        with a duty cycle between 0.25 and 0.5
%     available flags:
%      'CellType'  (should be a name, e.g. 'LP')
%      'Category' (0, 1, 2, 3)
%     the rest are ranges (e.g. [0.9, 1.0])
%      'BurstFreq', 'AutoCorr', 'SpikesPerBurst', ...
%      'BurstSpikeFreq', 'DutyCycle', 'SlowWaveAmp'
%  OUTPUTS:
%    MatchList: an array of structures
%       MatchList(n).matchData: structure with fields that record
%                               values that were searched for

% MatchList = GetSpecifiedAnalysis(Experiments, varargin)
%  INPUT ARGUMENTS:
%     -Experiments:  array of structures containing experiment data
%     -varargin:  specify properties to narrow down results.  All
%            are individually optional, but must specify at least one.
%        Condition:  a string denoting cell-type (e.g. 'GM')
%        Conductances: a two-element array, [g_syn, g_h]
%        Category:  an integer number, (0-3) corresponding to network type
%        AutoCorr:  a float number greater than zero and < 1
%  OUTPUT:
%     -MatchList:  array of structures containing matching experiment data

[cellType, category, specList, searchList] = getSpecifications(varargin{:});

MatchList = [];

flagList = {'CellType', 'Category', 'BurstFreq', 'AutoCorr', ...
	    'SpikesPerBurst', 'BurstSpikeFreq', 'DutyCycle', 'SlowWaveAmp'};

for n=1:length(Experiments)
  anList = Experiments(n).Analysis;
  for m=1:length(anList)
    an = anList(m);
    if(length(cellType) > 0 && ~strcmp(stripNums(an.Condition), ...
				       lower(cellType)))
      continue
    end
    if(isfinite(category) && an.Cat ~= category)
      continue
    end
    
    matchData = []; ind = 3;
    
    val = an.CellReal.Burst.Freq; ind = 3;
    if ~InRange(val, specList{ind})
      continue
    elseif(searchList{ind})
      matchData.(flagList{ind}) = val;
    end
    val = an.CellReal.SlowWave.Corr; ind = 4;
    if ~InRange(val, specList{ind})
      continue
    elseif(searchList{ind})
      matchData.(flagList{ind}) = val;
    end
    val = an.CellReal.Burst.SpikesPerBurst.Mean; ind = 5;
    if ~InRange(val, specList{ind})
      continue
    elseif(searchList{ind})
      matchData.(flagList{ind}) = val;
    end
    val = an.CellReal.Burst.SpikeFreq; ind = 6;
    if ~InRange(val, specList{ind})
      continue
    elseif(searchList{ind})
      matchData.(flagList{ind}) = val;
    end
    val = an.CellReal.Burst.DutyCycle; ind = 7;
    if ~InRange(val, specList{ind})
      continue
    elseif(searchList{ind})
      matchData.(flagList{ind}) = val;
    end
    val = mean(an.CellReal.SlowWave.Amplitudes); ind = 8;
    if ~InRange(val, specList{ind})
      continue
    elseif(searchList{ind})
      matchData.(flagList{ind}) = val;
    end
    
    an.FolderNum = Experiments(n).FolderNum;
    an.ID = sprintf('%u_%02u_%s', an.FolderNum, an.ExpNum, an.Condition);

    an.matchData = matchData;
    MatchList = [MatchList, an];
  end
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = getSpecifications(varargin)
if(length(varargin) < 2 || mod(length(varargin), 2) == 1)
  error('Invalid number of input arguments')
end

flagList = {'CellType', 'Category', 'BurstFreq', 'AutoCorr', ...
	    'SpikesPerBurst', 'BurstSpikeFreq', 'DutyCycle', 'SlowWaveAmp'};
flagList = lower(flagList);
specList = {'', NaN, [0, Inf], [-Inf, Inf], ...
	    [0, Inf], [0, Inf], [0, Inf], [0, Inf]};
searchList = {false, false, false, false, ...
	      false, false, false, false};
cellType = '';
category = NaN;

n = 1;
while(n < length(varargin))
  flag = varargin{n};
  if ~ischar(flag)
    error(['Error with argument %d: \n', ...
	   'Search specifications must be strings.'], n + 1)
  end
  lowFlag = lower(flag);
  flagInd = find(strcmp(flagList, lowFlag));
  if(length(flagInd) == 0)
    error('Unknown search specification: %s', lowFlag)
  end
  n = n + 1;
  if(flagInd == 1)
    cellType = varargin{n};
    if ~ischar(flag)
      error('CellType must be a string.')
    elseif ~(strcmp(lower(cellType), 'gm') || ...
	     strcmp(lower(cellType), 'dg') || ...
	     strcmp(lower(cellType), 'lp') || ...
	     strcmp(lower(cellType), 'pd'))
      error('Unknown cell type: %s', cellType)
    end
  elseif(flagInd == 2)
    category = varargin{n};
  else
    specList{flagInd} = varargin{n};
    if(length(specList{flagInd}) ~= 2)
      error('Search specifications must be a range.')
    end
    searchList{flagInd} = true;
  end
  n = n + 1;
end
varargout = {cellType, category, specList, searchList};
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stripStr = stripNums(conditionStr)
stripStr = regexp(conditionStr, '[a-zA-Z]*', 'match');
if(length(stripStr) == 0)
  error('Cell type %s is invalid:  must contain letters', conditionStr)
else
  stripStr = lower(stripStr{1});
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isInRange = InRange(x, range)
if( x < range(1) || x > range(2) )
  isInRange = false;
else
  isInRange = true;
end
return