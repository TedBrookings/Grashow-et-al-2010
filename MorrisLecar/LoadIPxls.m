function Intrinsics = LoadIPxls(xlsFileName)
if(nargin < 1)
  if(ispc)
    xlsFileName = ['C:\Documents and Settings\Rachel\My Documents\' ...
		   'spectra\ML analysis\Ted_ML_IPs.xls'];
  else
    xlsFileName = '/mnt/dwidget/ML analysis/Ted_ML_IPs.xls';
    %xlsFileName = '/mnt/dwidget/compen_ML_IPs.xls'
  end
end

%Turn off stupid "warning"
warning off MATLAB:xlsread:Mode
%Read in data
[~, ~, Raw] = xlsread(xlsFileName, 'Main', '', 'basic');


%Loop through data rows and pack info into structure
Intrinsics = [];
for Row = 2:size(Raw, 1)
  ID = getVal(Raw, Row, 'A');
  if(isempty(ID) || (length(ID) == 1 && isnan(ID)))
    continue
  end
  Use = getVal(Raw, Row, 'N');  
  if(~strcmpi(Use, 'Y'))
    continue
  end  
  
  IP.ID = ID;
  ParseInds = strfind(ID, '_');
  IP.FolderNum = ID(1:(ParseInds(1)-1));
  IP.ExpNum = ID((ParseInds(1)+1):(ParseInds(2)-1));
  IP.Condition = ID((ParseInds(2)+1):end);
  IP.BaseCond = StripNums(IP.Condition);
  
  IP.VRest = getVal(Raw, Row, 'C');
  IP.VRestChange = getVal(Raw, Row, 'D');
    
  IP.R = getVal(Raw, Row, 'F');
  IP.RChange = getVal(Raw, Row, 'G');
  
  IP.VThresh = getVal(Raw, Row, 'I');
  IP.VThreshChange = getVal(Raw, Row, 'J');
  
  IP.FISlope = getVal(Raw, Row, 'L');
  IP.FISlopeChange = getVal(Raw, Row, 'M');
  
  IP.SpikeRate1nA = getVal(Raw, Row, 'R');
  
  IP.SpikeHeight = getVal(Raw, Row, 'S');
  
  Intrinsics = [Intrinsics, IP];
end
return

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
if(isempty(StripStr))
  ErrStr = fprintf('Experiment type %s is invalid:  %s', ConditionStr, ...
		   'must contain letters');
  error(ErrStr)
else
  StripStr = StripStr{1};
end
return