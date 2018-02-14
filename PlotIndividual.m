function PlotIndividual(DirName)
%close all;
Records = OpenParamFile(DirName);
load([DirName, '/Analysis.mat']);  %load Analysis

%First find all the different categories ('ptx', '5ht', etc)
[CatList, IsMorrisLecar] = GetCatList(Records);
if(IsMorrisLecar == -1)
  PlotGMGM(DirName, Analysis, Records, CatList)
else
  PlotML(DirName, Analysis, Records, CatList)
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotGMGM(DirName, Analysis, Records, CatList)
ValueString = 'HalfCenter.Freq';  %Tells which values to map
PlotTrendBucking = true;

NumCat = length(CatList);
%plot 3d maps for each category
ptxNum = find(strcmp(CatList, 'ptx'));
CatRecords = GetCatFiles('ptx', Records);
ptxMap = MakeMap(DirName, 'ptx', CatRecords, Analysis, ValueString);
MapList = [];
for n = 1:NumCat
  disp(sprintf('Drawing map for %s (%d of %d)', CatList{n}, n, NumCat))
  CatRecords = GetCatFiles(CatList{n}, Records);
  
  if(n == ptxNum)
    Map = ptxMap;
    DrawMap(Map, ValueString);
  else
    Map = MakeMap(DirName, CatList{n}, CatRecords, Analysis, ValueString);
    DrawMap(Map, ValueString);
    DrawMap(Map, ValueString, ptxMap);
  end
  MapList = [MapList, Map];
  
  %Note the 'ptx' (control) map for later
  if(strcmp(CatList{n}, 'ptx'))
    ControlMapNum = n;
  end
end
clear ptxMap;

%Next plot the networks that slowed under each modulator type
NonControl = [1:(ControlMapNum-1), (ControlMapNum+1):NumCat];
ControlMap = MapList(ControlMapNum);
for n = NonControl
  disp(sprintf('Drawing decrease for %s (%d of %d)', ...
	       CatList{n}, n, NumCat-1))
  ModMap = MapList(n);
  if(PlotTrendBucking)
    DrawDecrease(ControlMap, ModMap, DirName);
    DrawDecrease(ControlMap, ModMap, DirName, 'ptx');
    DrawDecrease(ControlMap, ModMap, DirName, 'mod');
  else
    DrawDecrease(ControlMap, ModMap, DirName);
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotML(DirName, Analysis, Records, CatList)
ValueString = 'CellReal.Burst.Freq';  %Tells which values to map
PlotTrendBucking = true;

NumCat = length(CatList);
%plot 3d maps for each category
MapList = [];
for n = 1:NumCat
  disp(sprintf('Drawing map for %s (%d of %d)', CatList{n}, n, NumCat))
  CatRecords = GetCatFiles(CatList{n}, Records);
  
  Map = MakeMap(DirName, CatList{n}, CatRecords, Analysis, ValueString);
  DrawMap(Map, ValueString);
  DrawSpectrumMap(Map, 'cell');
  DrawSpectrumMap(Map, 'model');
  MapList = [MapList, Map];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Records = OpenParamFile(DirName)
FileList = dir(DirName);
FileList = GetTextFiles(FileList);

if(length(FileList) ~= 1)
  error(['Error in OpenParamFile.  Must specify directory with' ...
	 ' exactly one .txt file.']);
end

FileName = [DirName, '/', FileList.name];
fid = fopen(FileName);
if(fid < 0)
  error(['Error opening ', FileName]);
end

dummy = fgets(fid);  %contains the header description

Records = [];

while(1)
  Temp = fscanf(fid, '%d', 3);
  if(length(Temp) < 3)
    break;
  end
  Cat = fscanf(fid, '%s', 1);
  
  This.g_syn = Temp(1);
  This.g_h = Temp(2);
  This.Cat = Cat;
  This.AbfNamePart = sprintf('%s_%.4d.abf', Cat, Temp(3));
  This.SmrNamePart = sprintf('%s_%.4d.smr', Cat, Temp(3));
  This.DatNamePart = sprintf('%s_%.4d.dat', Cat, Temp(3));
  This.MatNamePart = sprintf('%s_%.4d.mat', Cat, Temp(3));
  
  Records = [Records, This];
end

fclose(fid);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TextFileList = GetTextFiles(FileList)
TextFileList = [];
for n = 1:length(FileList)
  if(strfind(FileList(n).name, '.txt'))
    TextFileList = [TextFileList, FileList(n)];
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CatList, IsMorrisLecar] = GetCatList(Records)
CatList = {Records(1).Cat};
NumCats = 1;

IsMorrisLecar = 0;
for n = 2:length(Records)
  Ind = find(strcmp(CatList, Records(n).Cat));
  if(length(Ind) == 0)  %new category
    NumCats = NumCats + 1;
    CatList(NumCats) = {Records(n).Cat};
    
    IsMorrisLecar = IsMorrisLecarCat(CatList{NumCats}, IsMorrisLecar);
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CatRecords = GetCatFiles(Cat, Records)
CatRecords = [];
for n = 1:length(Records)
  if(strcmp(Records(n).Cat, Cat))
    CatRecords = [CatRecords, Records(n)];
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Map = MakeMap(DirName, Cat, CatRecords, Analysis, ValueString)
g_syn = [];
g_h = [];
Results = [];

for n = 1:length(CatRecords)
  MatchCheck = 0;
  AbfCatFileName = CatRecords(n).AbfNamePart;
  SmrCatFileName = CatRecords(n).SmrNamePart;
  DatCatFileName = CatRecords(n).DatNamePart;
  MatCatFileName = CatRecords(n).MatNamePart;
  g_syn = [g_syn, CatRecords(n).g_syn];
  g_h = [g_h, CatRecords(n).g_h];
  for m = 1:length(Analysis)
    if(length(strfind(Analysis(m).FileName, AbfCatFileName)) > 0)
      Results = [Results, Analysis(m)];
      MatchCheck = 1;
      break;
    elseif(length(strfind(Analysis(m).FileName, SmrCatFileName)) > 0)
      Results = [Results, Analysis(m)];
      MatchCheck = 1;
      break;
    elseif(length(strfind(Analysis(m).FileName, DatCatFileName)) > 0)
      Results = [Results, Analysis(m)];
      MatchCheck = 1;
      break;
    elseif(length(strfind(Analysis(m).FileName, MatCatFileName)) > 0)
      Results = [Results, Analysis(m)];
      MatchCheck = 1;
    end
  end
  if(MatchCheck == 0)
    ErrStr = sprintf(['Could not find %s (or .smr, .dat, .mat)\n', ...
		      'Mismatch between experiment list and data files'], ...
		     AbfCatFileName);
    error(ErrStr);
  end
end

Map.g_syn = g_syn;
Map.g_h = g_h;
Map.Results = Results;
Map.Type = Cat;

Ind = [];
Values = [];
Cats = [];
for n = 1:length(Map.Results)
  Values = [Values, GetMapValue(Map.Results(n), ValueString)];
  Cats = [Cats, Map.Results(n).Cat];
  if(MapFilter(Map.Results(n)))
    Ind = [Ind, n];
  end
end

Map.Values = Values;
Map.Cats = Cats;
Map.Ind = Ind;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IncludeInMap = MapFilter(Analysis)
IncludeInMap = 1;

if(Analysis.Cat ~= 3)
  IncludeInMap = 0;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MapValue = GetMapValue(Analysis, varargin)
if(nargin >= 2)
  ValueString = varargin{1};
  MapValue = eval(['Analysis.', ValueString]);
else
  %MapValue = Analysis.Cat;
  %MapValue = (Analysis.Cell0.Spike.Freq + Analysis.Cell1.Spike.Freq)/ 2;
  MapValue = Analysis.HalfCenter.Freq;
  %MapValue = Analysis.Cell0.Spike.Frequencies.CoefOfVar;
end
return
