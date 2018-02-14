function [netProps, netLabels] = GetNarrowNetworks()

netLabels = {'g_synNarrow', 'gh_Narrow'};

xlsData = getXlsData();

netProps = convertFormat(xlsData);
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
function xlsData = getXlsData()
if(ispc)
  xlsFileName = ['C:\Documents and Settings\Rachel\My Documents\' ...
		 'spectra\Narrow_selection_info.xls'];
else
  xlsFileName = '/mnt/dwidget/Narrow_selection_info.xls';
end

%Turn off stupid "warning"
warning off MATLAB:xlsread:Mode
%Read in data
[numFields, textFields, rawFields] = ...
    xlsread(xlsFileName, 'Sheet1', '', 'basic');
xlsData = [];
numLines = size(numFields, 1);
for n = 1:numLines
  temp.netID = textFields{n+1,1};
  ind = strfind(temp.netID, '_');
  temp.ID = temp.netID(1:(ind(3)-1));
  temp.folderNum = str2num(temp.netID(1:(ind(1)-1)));
  temp.expNum = str2num(temp.netID((ind(1)+1):(ind(2)-1)));
  temp.condition = temp.netID((ind(2)+1):(ind(3)-1));
  temp.cellType = StripNums(temp.condition);
  
  temp.g_syn = numFields(n,1);
  temp.g_h = numFields(n,2);
  
  xlsData = [xlsData, temp];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StripStr = StripNums(CellTypeStr)
StripStr = regexp(CellTypeStr, '[a-zA-Z]*', 'match');
if(length(StripStr) == 0)
  ErrStr = fprintf('Cell type %s is invalid:  %s', CellTypeStr, ...
		   'must contain letters');
  error(ErrStr)
else
  StripStr = StripStr{1};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function netProps = convertFormat(xlsData)
numNets = length(xlsData);
numFields = 2;
netProps = cell(1, numFields);
for n = 1:numFields
  tempList = [];
  for m = 1:numNets
    net = xlsData(m);
    temp.CellType = net.cellType;
    temp.ID = net.ID;
    if(n == 1)
      temp.Mean = net.g_syn;
    else
      temp.Mean = net.g_h;
    end
    tempList = [tempList, temp];
  end
  netProps{1,n} = tempList;
end
return