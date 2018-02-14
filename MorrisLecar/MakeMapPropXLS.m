function MakeMapPropXLS()
%options.netType = 1;   %Cell spiking
options.netType = 3;  %Bursting
%options.netType = 0.85;  %Highly-regular slow wave

options.numSynDivisions = 2;

xlsFile = '/mnt/dwidget/ML analysis/MapProperties.csv';

%% load the experiment data

experiments = getExperiments();
intrinsics = LoadIPxls;

%% Get lots of properties from the maps:
[mapProps, mapPropLabels] = getAllMapProps(experiments, options);

%% Make the spreadsheet table
table = makeSpreadSheet(mapProps, mapPropLabels, intrinsics);

%% write table to xlsFile
writeTable(table, xlsFile);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function experiments = getExperiments()
experimentsExists = evalin('base', 'exist(''ML_Experiments'')');
if(experimentsExists)
  experiments = evalin('base', 'ML_Experiments');
else
  experiments = LoadAllExperiments('MorrisLecarFolders.txt');
  assignin('base', 'ML_Experiments', experiments);
end

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
mapProps = [MapQuantities, {synSens}, {hSens}, LocalDensities];
mapPropLabels = [MapQuantitiesLabels, SensLabels, ...
		 DenseLabels];
if length(mapProps) ~= length(mapPropLabels)
  fprintf(2, 'Inconsistent number of map properties!\n')
  keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function table = makeSpreadSheet(mapProps, mapPropLabels, intrinsics)
IDList = {intrinsics.ID};
numCells = length(IDList);
numProps = length(mapProps);
mapPropIDList = {mapProps{1}.ID};

table = cell(numCells + 1, numProps + 1);
table{1,1} = 'cell ID';
for prop = 1:numProps
  col = prop + 1;
  label = mapPropLabels{prop};
  table{1, col} = label;
  
  mapPropID = {mapProps{prop}.ID};
  for cellNum = 1:numCells
    row = cellNum + 1;
    ID = IDList{cellNum};
    
    mapPropRow = find(strcmp(mapPropID, ID));
    if isempty(mapPropRow)
      value = 'NaN';
    elseif length(mapPropRow) == 1
      value = mapProps{prop}(mapPropRow).Mean;
    else
      error('Multiple data sets for %s %s', ID, ...
	    mapPropLabels{prop})
    end
    table{row, col} = value;
  end
end

for cellNum = 1:numCells
  row = cellNum + 1;
  ID = IDList{cellNum};
  table{row, 1} = ID;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeTable(table, csvFile)
%% write the spreadsheet to an xls file
%Turn off stupid "warning"
fid = fopen(csvFile, 'w');
[numRows, numCols] = size(table);
for row = 1:numRows
  for col = 1:numCols
    val = table{row,col};
    if col == numCols
      delim = '';
    else
      delim = ',';
    end
    
    if isempty(val)
      fprintf(fid, '%s', delim);
    elseif ischar(val)
      fprintf(fid, '%s%s', val, delim);
    else
      fprintf(fid, '%g%s', val, delim);
    end
  end
  fprintf(fid, '\r\n');
end
fclose(fid);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeTableDumb(table, xlsFile)
%% write the spreadsheet to an xls file
%Turn off stupid "warning"
warning off MATLAB:xlswrite:NoCOMServer
sheet = 'main';
[numRows, numCols] = size(table);
for row = 1:numRows
  for col = 1:numCols
    val = table{row,col};
    if isempty(val)
      continue
    end
    xlsCol = getXlsCol(col);
    range = sprintf('%s%d', xlsCol, row);
    [status, message] = xlswrite(xlsFile, val, sheet, range);
    if ~status
      fprintf(2, '%s\n', message.message)
      break
    end
  end
  if ~status
    break
  end
end

if status
  fprintf('%s successfully created.\n', xlsFile)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xlsCol = getXlsCol(col)
remainder = mod(col, 26);
xlsCol = char('A' + remainder - 1);
while col > 26
  col = round((col - remainder) / 26);
  remainder = mod(col, 26);
  xlsCol = [char('A' + remainder - 1), xlsCol];
end
return