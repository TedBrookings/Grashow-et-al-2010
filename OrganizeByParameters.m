function Organized = OrganizeByParameters(Experiments, OrgString)
if(nargin < 1 || length(Experiments) == 0)
  ExperimentsExists = evalin('base', 'exist(''Experiments'')');
  if(ExperimentsExists)
    Experiments = evalin('base', 'Experiments');
  else
    if(nargin < 2)
      OrgString = 'HalfCenter.Freq';
      ListFile = 'folder_names.txt';
    elseif(strcmp(OrgString, 'HalfCenter.Freq'))
      ListFile = 'folder_names.txt';
    else
      ListFile = 'MorrisLecarFolders.txt';
    end
      
    Experiments = LoadAllExperiments(ListFile);
    %Uncomment this line to save Experiments to the base workspace
    %assignin('base', 'Experiments', Experiments);
  end
elseif(nargin < 2)
  OrgString = 'HalfCenter.Freq';
end

OrgString = ['Freq = [Freq, Trial(m).', OrgString, '];'];

%Construct all the lists we want:
g_syn = [];
g_h = [];
Freq = [];
ExpNum = [];
Condition = {};
Analysis = [];

NumExperiments = length(Experiments);
for n = 1:NumExperiments
    Trial = Experiments(n).Analysis;
    NumTrials = length(Trial);
    for m = 1:NumTrials
        g_syn = [g_syn, Trial(m).g_syn];
        g_h = [g_h, Trial(m).g_h];
	eval(OrgString);  %Unless OrgString is altered, executes:
			  %   Freq = [Freq, Trial(m).HalfCenter.Freq];
        ExpNum = [ExpNum, Trial(m).ExpNum];
        %Temp = Trial(m).Condition;
        Condition = {Condition{:}, Trial(m).Condition};
	Analysis = [Analysis, Trial(m)];
    end
end

Organized = ReArrange(g_syn, g_h, Freq, ExpNum, Condition, Analysis);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Organized = ReArrange(g_syn, g_h, Freq, ExpNum, Condition, Analysis)
%First find all the condition types, and all the unique
%  (g_syn, g_h, ExpNum) combinations
ConditionList = {StripNums(Condition{1})};
g_syn_List = g_syn(1);
g_h_List = g_h(1);
ExpNum_List = ExpNum(1);
%Also, form lists of how the raw data matches up to the organized data
ConditionNums = zeros(size(Condition));
ConditionNums(1) = 1;
gNums = zeros(size(g_syn));
gNums(1) = 1;
for n = 2:length(Condition)
    Match = 0;
    Cond_n = StripNums(Condition{n});
    for m=1:length(ConditionList)
       if(strcmp(ConditionList{m}, Cond_n))
           Match = 1;
           ConditionNums(n) = m;
           break;
       end
    end
    if(Match == 0)
        ConditionList = {ConditionList{:}, Cond_n};
        ConditionNums(n) = length(ConditionList);
    end
    
    Match = 0;
    for m=1:length(g_syn_List)
      if(g_syn(n) == g_syn_List(m) & g_h(n) == g_h_List(m) ...
	 & ExpNum(n) == ExpNum_List(m))
        Match = 1;
        gNums(n) = m;
	break;
      end
    end
    if(Match == 0)
        g_syn_List = [g_syn_List, g_syn(n)];
        g_h_List = [g_h_List, g_h(n)];
        ExpNum_List = [ExpNum_List, ExpNum(n)];
        gNums(n) = length(g_syn_List);
    end
end

%Make a matrix to house every possible frequency
Num_Cond = length(ConditionList);
Num_g = length(g_syn_List);
FreqMat = repmat(NaN, Num_g, Num_Cond + 3);
FreqMat(:,1) = g_syn_List';
FreqMat(:,2) = g_h_List';
FreqMat(:,3) = ExpNum_List';

%Now add the raw data into the organized matrix
% and make a (temporary) matrix to house analysis structs
for n = 1:length(g_syn)
  Row = gNums(n);
  Col = ConditionNums(n);
  FreqMat(Row,Col+3) = Freq(n);
  AnalysisMat(Row, Col) = Analysis(n);
end

%Finally, put all this into the structure "Organized"
Organized.g_syn = g_syn_List;
Organized.g_h = g_h_List;
Organized.ExpNum = ExpNum_List;
for n = 1:length(ConditionList)
  TypeName = sprintf('Type_%s', ConditionList{n});
  Organized.(TypeName) = FreqMat(:,n+3);
  TypeAnalysisName = sprintf('Analysis_%s', ConditionList{n});
  Organized.(TypeAnalysisName) = AnalysisMat(:,n);
end
Organized.FreqMat = FreqMat;
Organized.Conditions = ConditionList;
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
