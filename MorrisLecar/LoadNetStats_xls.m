function NetStats = LoadNetStats_xls(xlsFileName)
WinDir = 'C:\Documents and Settings\Rachel\My Documents\spectra\ML analysis\';
NixDir = '/mnt/dwidget/ML analysis/';
if(ispc)
  DirName = WinDir;
else
  DirName = NixDir;
end

if(nargin < 1 || StringCheck(xlsFileName, 'LP'))
  xlsFileName = [DirName, 'LP_network_props.xls'];
  cellType = 'LP';
elseif(StringCheck(xlsFileName, 'PD'))
  xlsFileName = [DirName, 'PD_network_props.xls'];
  cellType = 'PD';
else
  error('Invalid file requested.');
end
%Turn off stupid "warning"
warning off MATLAB:xlsread:Mode
%Read in data
[Num, Text, Raw] = xlsread(xlsFileName, 1, '', 'basic');
numRows = size(Raw, 1);

%Loop through data rows and pack info into structure
NetStats = [];
for Row = 2:numRows
  ID = getVal(Raw, Row, 'A');
  if(length(ID) == 0 || (length(ID) == 1 && isnan(ID)))
    continue
  end
  NS.ID = ID;
  ParseInds = strfind(ID, '_');
  NS.FolderNum = ID(1:(ParseInds(1)-1));
  NS.ExpNum = ID((ParseInds(1)+1):(ParseInds(2)-1));
  NS.Condition = ID((ParseInds(2)+1):end);
  NS.BaseCond = StripNums(NS.Condition);
  
  NS.CyclePeriod = getVal(Raw, Row, 'H');  %s
  NS.BurstFrequency = getVal(Raw, Row, 'I');  %Hz
  NS.Duration = getVal(Raw, Row, 'J');
  NS.DutyCycle = getVal(Raw, Row, 'K');
  NS.SpikesPerBurst = getVal(Raw, Row, 'M');
  NS.SpikeFreq = getVal(Raw, Row, 'N');

  NS.LPOnPhase = getVal(Raw, Row, 'AH');
  NS.LPOffPhase = getVal(Raw, Row, 'AI');
  NS.PYOnPhase = getVal(Raw, Row, 'AV');
  NS.PYOffPhase = getVal(Raw, Row, 'AW');
  NS.PDOffPhase = getVal(Raw, Row, 'BD');
  
  if strcmp(cellType, 'LP')
    NS.Phase = NS.LPOnPhase;
  else
    NS.Phase = NS.PDOffPhase;
  end
  
  NetStats = [NetStats, NS];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = getVal(table, rowNum, colStr)
colNum = 0;
for n = 1:length(colStr)
  val = 1 + (colStr(n) - 'A');
  colNum = colNum + val * 26^(length(colStr) - n);
end
val = table{rowNum, colNum};
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