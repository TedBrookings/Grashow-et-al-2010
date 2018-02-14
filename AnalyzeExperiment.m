function varargout = AnalyzeExperiment(DirName)
tic
if(ispc)
  Slash = '\';
  AltSlash = '/';
else
  Slash = '/';
  AltSlash = '\';
end
DefaultDir = 'C:\Documents and Settings\Rachel\My Documents\spectra\';
LinuxDir = '/mnt/dwidget/';

%Add trailing Slash
if(DirName(end) ~= Slash)
  DirName = [DirName, Slash];
end
[Records, IsMorrisLecar] = OpenParamFile(DirName);

DirStruct = dir(DirName);

FileList = [];
OrderInd = zeros(size(Records));
for n = 1:length(DirStruct)
  %Check if the file in DirStruct(n).name is in the list of records
  %we want to analyze:
  RecordInd = IsInRecordsList(DirStruct(n).name, Records);
  if(RecordInd > 0)
    if(Records(RecordInd).Matched == true)  %make sure that it can only
      continue;                             %be matched once (in case
    else                                    %there is a .abf AND a .smr
      Records(RecordInd).Matched = true;
    end
    
    FileStruct = DirStruct(n);
    FileStruct.g_syn = Records(RecordInd).g_syn;
    FileStruct.g_h = Records(RecordInd).g_h;
    FileStruct.Condition = Records(RecordInd).Cat;
    FileList = [FileList, FileStruct];
    OrderInd(RecordInd) = length(FileList);
  end
end
MissedInd = find(OrderInd == 0);
if(length(MissedInd) > 0)
  fprintf(2, 'Missed file number(s) in %s:', DirName);
  for Ind = MissedInd
    fprintf(2, ' (%s %g)', Records(Ind).Cat, Records(Ind).FileNum)
  end
  fprintf(2, '%s', '\n');
  error(['Error analyzing ', DirName])
end
FileList = FileList(OrderInd);  %Analyze in order of records list

NumFiles = length(FileList);
if(NumFiles == 0)
  error(['No files to analyze in ', DirName]);
end

OutFileName = [DirName, 'Analysis.mat'];
ContinueFileName = [DirName, 'Analysis_Continue.mat'];
[Analysis, Experiments, Start_n] = CheckContinue(ContinueFileName);
DoContinue = true;
n = Start_n;
%One would wish this worked, but it doesn't:
%CleanObj = onCleanup(@()SaveResults(OutFileName, Analysis, ...
%				    Experiments, n, DoContinue));

for n = Start_n:NumFiles
%  CleanObj = onCleanup(@()SaveResults(OutFileName, Analysis, ...
%				      Experiments, n, DoContinue));
  try
    Experiments = [Experiments, GetExpNum(FileList(n).name)];
    FileName = [DirName, FileList(n).name];
    fprintf('Analyzing %s (%d of %d)\n', FileName, n, NumFiles)
    if(IsMorrisLecar)
      OutStruct = AnalyzeML(FileName);
    else
      OutStruct = RunAnalyze(FileName, 0);
    end
    if(strcmp(DirName, LinuxDir))
      OutStruct.FileName = [DefaultDir, FileList(n).name];
    else
      OutStruct.FileName = FileName;
    end
    OutStruct.ExpNum = Experiments(end);
    OutStruct.g_syn = FileList(n).g_syn;
    OutStruct.g_h = FileList(n).g_h;
    OutStruct.Condition = FileList(n).Condition;
    Analysis = [Analysis, OutStruct];
  catch
    ErrStruct = lasterror;
    if(n > Start_n)
      SaveResults(OutFileName, ContinueFileName, DoContinue, ...
		  Analysis, Experiments, n);
    end
    rethrow(ErrStruct)
  end
end
DoContinue = false;
%CleanObj = onCleanup(@()SaveResults(OutFileName, Analysis, ...
%				    Experiments, n, DoContinue));  
ExpNum = Experiments(1);
if(sum(Experiments ~= ExpNum) > 0)
  ExpNum = -1;
end

if(nargout <= 1)
  varargout = {Analysis};
elseif(nargout == 2)
  varargout = {Analysis, ExpNum};
end

SaveResults(OutFileName, ContinueFileName, DoContinue, ...
	    Analysis, Experiments, n);
Toc
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Records, IsMorrisLecar] = OpenParamFile(DirName)

FileList = dir(DirName);
FileList = GetTextFiles(FileList);

if(length(FileList) ~= 1)
  error(['Error in OpenParamFile.  Must specify directory with' ...
	 ' exactly one .txt file.']);
end

FileName = [DirName, FileList.name];

fid = fopen(FileName);
if(fid < 0)
  error(['Error opening ', FileName]);
end

dummy = fgets(fid);  %contains the header description

Records = [];
IsMorrisLecar = 0;

while(1)
  Temp = fscanf(fid, '%d', 3);
  if(length(Temp) < 3)
    break;
  end
  Cat = fscanf(fid, '%s', 1);
  
  This.g_syn = Temp(1);
  This.g_h = Temp(2);
  This.Cat = Cat;
  This.FileNum = Temp(3);
  
  IsMorrisLecar = IsMorrisLecarCat(Cat, IsMorrisLecar);
  
  This.AbfNamePart = sprintf('%s_%.4d.abf', Cat, Temp(3));
  This.SmrNamePart = sprintf('%s_%.4d.smr', Cat, Temp(3));
  This.DatNamePart = sprintf('%s_%.4d.dat', Cat, Temp(3));
  This.MatNamePart = sprintf('%s_%.4d.mat', Cat, Temp(3));
  This.Matched = false;
  
  Records = [Records, This];
end
IsMorrisLecar = (IsMorrisLecar > 0);

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
function RecordInd = IsInRecordsList(FileName, Records)
RecordInd = 0;
for n = 1:length(Records)
  if(length(strfind(FileName, Records(n).AbfNamePart)) > 0)
    RecordInd = n;
    return
  elseif(length(strfind(FileName, Records(n).SmrNamePart)) > 0)
    RecordInd = n;
    return
  elseif(length(strfind(FileName, Records(n).DatNamePart)) > 0)
    RecordInd = n;
    return
  elseif(length(strfind(FileName, Records(n).MatNamePart)) > 0)
    RecordInd = n;
    return
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExpNum = GetExpNum(FileName)
Ind = strfind(FileName, '_');
if(length(Ind) >= 2);
  Range = (Ind(1)+1):(Ind(2)-1);
  ExpNum = str2double(FileName(Range));
else
  ExpNum = -1;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Analysis, Experiments, Start_n] = CheckContinue(ContinueFileName)
DoContinue = false;
try
  load(ContinueFileName);
catch
  DoContinue = false;
end

if(DoContinue)
  disp(sprintf('Resuming from file# %g', Start_n))
else
  Analysis = [];
  Experiments = [];
  Start_n = 1;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveResults(OutFileName, ContinueFileName, DoContinue, ...
		     Analysis, Experiments, n)
if(DoContinue)
  if(n > 1)  %Try to save any actual progress
    Start_n = n;
    Analysis = Analysis(1:(n - 1));
    Experiments = Experiments(1:(n - 1));
    DoContinue = true;
    try
      save(ContinueFileName, 'Analysis', 'Experiments', 'Start_n', ...
	   'DoContinue');
    catch
      disp(sprintf('Warning, could''t save %s', ContinueFileName))
    end
  end
else
  ExpNum = Experiments(1);
  if(sum(Experiments ~= ExpNum) > 0)
    ExpNum = -1;
  end
  try
    save(OutFileName, 'Analysis', 'ExpNum');
    Saved = true;
  catch
    disp(sprintf('Warning, could''t save %s', OutFileName))
    Saved = false;
  end
  if(Saved)  %delete the continue file, if it exists
    fid = fopen(ContinueFileName, 'r');
    if(fid > 0)
      fclose(fid);
      delete(ContinueFileName);
    end
  end
end

return