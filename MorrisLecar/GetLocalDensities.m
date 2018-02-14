function [Densities, DenseNames] = GetLocalDensities(Experiments, ...
						  NetType, NumDivide)

if(nargin < 3)
  NumDivide = 2;
end
CellTypeList = GetCellTypeList(Experiments);

UseHardCutoff = false;
DivideTwo = false;

if(abs(NetType - round(NetType)) < 1e-9 || UseHardCutoff)
  FillStyle = 1;
else
  FillStyle = 2;
end


if(NumDivide == 2)
  DenseNames = {'LL', 'HL', ...
		'LH', 'HH'};
  SynSets = {[10, 25, 40, 55], [55, 70, 85, 100], ...
	     [10, 25, 40, 55], [55, 70, 85, 100]};
  SynCoefs = {[1.0, 1.0, 1.0, 0.5], [0.5, 1.0, 1.0, 1.0], ...
	      [1.0, 1.0, 1.0, 0.5], [0.5, 1.0, 1.0, 1.0]};
  HSets = {[10, 25, 40, 55], [10, 25, 40, 55], ...
	   [55, 70, 85, 100], [55, 70, 85, 100]};
  HCoefs = {[1.0, 1.0, 1.0, 0.5], [1.0, 1.0, 1.0, 0.5], ...
            [0.5, 1.0, 1.0, 0.5], [0.5, 1.0, 1.0, 1.0]};
  Densities = {[], [], ...
	       [], []};
  IDList = {};
else
  DenseNames = {'LL', 'ML', 'HL', ...
		'LM', 'MM', 'HM', ...
		'LH', 'MH', 'HH'};
  SynSets = {[10, 25, 40], [40, 55, 70], [70, 85, 100], ...
	     [10, 25, 40], [40, 55, 70], [70, 85, 100], ...
	     [10, 25, 40], [40, 55, 70], [70, 85, 100]};
  SynCoefs = {[1.0, 1.0, 1.0/3], [2.0/3, 1.0, 2.0/3], [1.0/3, 1.0, 1.0], ...
	      [1.0, 1.0, 1.0/3], [2.0/3, 1.0, 2.0/3], [1.0/3, 1.0, 1.0], ...
	      [1.0, 1.0, 1.0/3], [2.0/3, 1.0, 2.0/3], [1.0/3, 1.0, 1.0]};
  HSets = {[10, 25, 40], [10, 25, 40], [10, 25, 40], ...
	   [40, 55, 70], [40, 55, 70], [40, 55, 70], ...
	   [70, 85, 100], [70, 85, 100], [70, 85, 100]};

  HCoefs = {[1.0, 1.0, 1.0/3], [1.0, 1.0, 1.0/3], [1.0, 1.0, 1.0/3], ...
	    [2.0/3,1.0,2.0/3], [2.0/3,1.0,2.0/3], [2.0/3,1.0,2.0/3], ...
	    [1.0/3, 1.0, 1.0], [1.0/3, 1.0, 1.0], [1.0/3, 1.0, 1.0]};
  Densities = {[], [], [], ...
	       [], [], [], ...
	       [], [], []};
  IDList = {};
end

if(NetType == round(NetType))
  Flag = 'Category';
else
  Flag = 'AutoCorr';
  NetType = [NetType, Inf];
end

for n=1:length(CellTypeList)
  CellType = CellTypeList{n};
  
  if(FillStyle == 1)
    %Get the networks of type CellType, with network properties NetType
    NetList = GetSpecifiedAnalysis(Experiments, 'CellType', CellType, ...
				   Flag, NetType);
  else
    %Get all networks of type CellType
    NetList = GetSpecifiedAnalysis(Experiments, 'CellType', CellType);
  end
  NumCells = length(unique({NetList.ID}));
  if(NumCells < 4)
    continue
  end

  %Loop through all the networks that meet our criterion
  NumNetworks = length(NetList);
  for m = 1:NumNetworks
    Network_m = NetList(m);
    g_syn_m = Network_m.g_syn;
    g_h_m = Network_m.g_h;
    ID_m = Network_m.ID;
    [MatchID, IDInd] = ismember(ID_m, IDList);
    if(~MatchID)
      %Execute this block if this is the first network for a given cell
      %Make an empty structure to hold density info
      IDList = {IDList{:}, ID_m};
      IDInd = length(IDList);
      Temp.CellType = CellType;
      Temp.ID = ID_m;
      Temp.Mean = 0;
      %Append a copy of that structure to the end of each Densities location
      for k = 1:length(Densities)
	Densities{k} = [Densities{k}, Temp];
      end
    end
    if(FillStyle == 1)
      Val_m = 1.0;
    else
      Val_m = Network_m.CellReal.SlowWave.Corr;
    end
    for k = 1:length(Densities)
      %For each Densities location, see if this network has any overlap
      [SynMatch, SynInd] = ismember(g_syn_m, SynSets{k});
      [HMatch, HInd] = ismember(g_h_m, HSets{k});
      if(SynMatch && HMatch)
        %If so, calculate the coefficient, and add it to the density
        Coef = SynCoefs{k}(SynInd) * HCoefs{k}(HInd);
        Densities{k}(IDInd).Mean = Densities{k}(IDInd).Mean ...
	    + Coef * Val_m;
        %Used to break here, but now a point may overlap more than one location
        %break
      end
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CellTypeList = GetCellTypeList(Experiments)
CellTypeList = {};
for n = 1:length(Experiments)
  Temp = Experiments(n).ConditionList;
  for m = 1:length(Temp)
    Temp{m} = StripNums(Temp{m});
  end
  CellTypeList = {CellTypeList{:}, Temp{:}};
end
CellTypeList = unique(CellTypeList);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StripStr = StripNums(ConditionStr)
StripStr = regexp(ConditionStr, '[a-zA-Z]*', 'match');
if(length(StripStr) == 0)
  ErrStr = fprintf('Experiment type %s is invalid:  %s', ConditionStr, ...
		   'must contain letters');
  error(ErrStr)
else
  StripStr = StripStr{1};
end
return