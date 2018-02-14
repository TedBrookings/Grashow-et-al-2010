%DEPRICATED!!!
function Organized = OrganizeML(Experiments)
%DEPRICATED!!!
if(nargin < 1 || length(Experiments) == 0)
  ExperimentsExists = evalin('base', 'exist(''Experiments'')');
  if(ExperimentsExists)
    Experiments = evalin('base', 'Experiments');
  else
    ListFile = 'MorrisLecarFolders.txt';
    Experiments = LoadAllExperiments(ListFile);
    %Uncomment this line to save Experiments to the base workspace
    %assignin('base', 'Experiments', Experiments);
  end
elseif(nargin < 2)
  OrgString = 'HalfCenter.Freq';
end

ConditionList = GetConditionList(Experiments);
Organized = [];
for n = 1:length(ConditionList)
  TempS.Condition = ConditionList{n};
  TempS.g_syn = [];
  TempS.g_h = [];
  TempS.Analysis = {};
  Organized = [Organized, TempS];
end

Organized = FillOrganizedStruct(Experiments, ConditionList, Organized);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ConditionList = GetConditionList(Experiments)
ConditionList = {};
for n = 1:length(Experiments)
  Temp = Experiments(n).ConditionList;
  for m = 1:length(Temp)
    Temp{m} = StripNums(Temp{m});
  end
  ConditionList = {ConditionList{:}, Temp{:}};
end
ConditionList = unique(ConditionList);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Organized = FillOrganizedStruct(Experiments, ConditionList, ...
					 Organized)
for n = 1:length(Experiments)
  Exp = Experiments(n);
  for m = 1:length(Exp.Analysis)
    An = Exp.Analysis(m);
    Cond = StripNums(An.Condition);
    OInd = find(strcmp(ConditionList, Cond));
    Org = Organized(OInd);
    gInd = find(Org.g_syn == An.g_syn & Org.g_h == An.g_h);
    if(length(gInd) == 0)
      Organized(OInd).g_syn = [Org.g_syn, An.g_syn];
      Organized(OInd).g_h = [Org.g_h, An.g_h];
      Organized(OInd).Analysis = {Org.Analysis{:}, []};
      gInd = length(Org.g_syn) + 1;
    end

    Organized(OInd).Analysis{gInd} = [Organized(OInd).Analysis{gInd}, An];
  end
end
  
return