function varargout = RunIntrinsic
GetSubDir = true;
if(ispc)
  BaseDir = 'C:\My Documents'
  Slash = '\';
else
  %BaseDir = '/home/ted/data/rachel';
  BaseDir = '/mnt/dwidget2';
  Slash = '/';
end
IVDir = [BaseDir, Slash, 'IV curves', Slash];
RampDir = [BaseDir, Slash, 'Ramps', Slash];
FIDir = [BaseDir, Slash, 'FI curves', Slash];

SeparateTop_Bot = true;

IVExists = evalin('base', 'exist(''IVList'')');
if(IVExists)
  IVList = evalin('base', 'IVList');
else
  JointFiles = false;
  IVList = GetExpList(IVDir, GetSubDir, SeparateTop_Bot, JointFiles);
  IVList = AnalyzeList(IVList, 'IV', @GetIVCurve);
  assignin('base', 'IVList', IVList);
end

RampExists = evalin('base', 'exist(''RampList'')');
if(RampExists)
  RampList = evalin('base', 'RampList');
else
  JointFiles = true;
  RampList = GetExpList(RampDir, GetSubDir, SeparateTop_Bot, JointFiles);
  RampList = AnalyzeList(RampList, 'Ramp', @GetRampData);
  assignin('base', 'RampList', RampList);
end

FIExists = evalin('base', 'exist(''FIList'')');
if(FIExists)
  FIList = evalin('base', 'FIList');
else
  JointFiles = true;
  FIList = GetExpList(FIDir, GetSubDir, SeparateTop_Bot, JointFiles);
  FIList = AnalyzeList(FIList, 'FI', @GetFICurve);
  assignin('base', 'FIList', FIList);
end

MakeIntrinsicPlots(IVList, RampList, FIList);

switch(nargout)
 case 0, varargout = {};
 case 1, varargout = {IVList};
 case 2, varargout = {IVList, RampList};
 case 3, varargout = {IVList, RampList, FIList};
 otherwise, error('Too many output arguments');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExpList = GetExpList(Dir, GetSubDir, SeparateTop_Bot, JointFiles)
FileList = GetFileList(Dir, GetSubDir, SeparateTop_Bot, JointFiles);

StructNames = unique({FileList.StructName});
BlankExp.ID = '';
BlankExp.Book = 0;
BlankExp.ExpNum = 0;
if(SeparateTop_Bot)
  BlankExp.WhichCell = '';
end
for n = 1:length(StructNames)
  BlankExp.Files.(StructNames{n}) = '';
end

RawIDList = {FileList.ID};
IDList = {};
ExpList = [];
for n = 1:length(RawIDList)
  This = RawIDList{n};
  [IsInList, Ind] = ismember(This, IDList);
  if(IsInList)
    ThisFile = FileList(n);
    ExpList(Ind).Files.(ThisFile.StructName) = ThisFile.FileName;
  else
    IDList = {IDList{:}, This};
    ThisExp = BlankExp;
    ThisFile = FileList(n);
    ThisExp.Files.(ThisFile.StructName) = ThisFile.FileName;
    ThisExp.Book = ThisFile.Book;
    ThisExp.ExpNum = ThisFile.ExpNum;
    if(SeparateTop_Bot)
      ThisExp.WhichCell = ThisFile.WhichCell;
    end
    ThisExp.ID = This;
    ExpList = [ExpList, ThisExp];
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FileList = GetFileList(Dir, GetSubDir, SeparateTop_Bot, ...
				JointFiles)
DirList = dir(Dir);
if(length(DirList)==0)
  error(['Directory ', Dir, ' does not exist.'])
end
if(strcmp(DirList(1).name, '.'))
  Start_n = 3;
else
  Start_n = 1;
end

if(ispc)
  Slash = '\';
else
  Slash = '/';
end

FileList = [];
NumDirList = length(DirList);
n = Start_n;
for n=Start_n:NumDirList
  if(DirList(n).isdir)
    if(GetSubDir)
      NewDir = [Dir, DirList(n).name, Slash];
      FileList = [FileList, GetFileList(NewDir, GetSubDir, ...
					SeparateTop_Bot, JointFiles)];
    end
    continue
  end
  Name = DirList(n).name;
  if(~JointFiles)
    ThisFile = GetFileInfo(Dir, Name, SeparateTop_Bot);
  else
    ThisFile = GetFileInfo(Dir, Name, SeparateTop_Bot, 'top');
    FileList = [FileList, ThisFile];
    ThisFile = GetFileInfo(Dir, Name, SeparateTop_Bot, 'bot');
  end     
  FileList = [FileList, ThisFile];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ThisFile = GetFileInfo(Dir, Name, SeparateTop_Bot, WhichCell)
		  
FileName = [Dir, Name];
Ind = sort([strfind(Name, '_'), strfind(Name, '.')]);

Book = sscanf(Name, '%g', 1);
ExpNum = sscanf(Name(Ind(1):end), '%g', 1);
Book_Exp = Name(1:(Ind(2)-1));
Condition = Name((Ind(2)+1):(Ind(3)-1));
if(nargin == 3)
  WhichCell = Name((Ind(3)+1):(Ind(3)+3));
end

CrazyNum = (Book_Exp(end) == 'a');

if(SeparateTop_Bot)
  ID = [Book_Exp, '_', WhichCell];
  StructName = sprintf('file_%s', Condition);
else
  ID = Book_Exp;
  StructName = sprintf('file_%s_%s', WhichCell, Condition);
end
  
ThisFile.FileName = FileName;
ThisFile.Book = Book;
ThisFile.ExpNum = ExpNum;
ThisFile.CrazyNum = CrazyNum;
ThisFile.Condition = Condition;
ThisFile.WhichCell = WhichCell;
ThisFile.Book_Exp = Book_Exp;
ThisFile.StructName = StructName;
ThisFile.ID = ID;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutList = AnalyzeList(ExpList, ListType, ListFunc)
OutList = [];
FNames = fieldnames(ExpList(1).Files);
NumExp = length(ExpList);
tic
for n = 1:NumExp
  This = ExpList(n);
  for m = 1:length(FNames)
    FileName = This.Files.(FNames{m});
    try
      if(length(FileName) == 0)
	This.(ListType).(FNames{m}) = ListFunc(); %returns a blank
	continue;
      end
      %disp(FileName)
      [t, v, I, Okay] = GetCleanData(FileName, This);
      if(~Okay)
	This.Files.(FNames{m}) = '';
	This.(ListType).(FNames{m}) = ListFunc(); %returns a bunch of NaNs
      else
	This.(ListType).(FNames{m}) = ListFunc(t, v, I);
      end
    catch
      disp(['Error analyzing ', FileName, ' ', This.WhichCell])
      rethrow(lasterror)
    end
  end
  OutList = [OutList, This];
  t = toc;
  ProgStr = sprintf('%g of %g %s files completed, %g min remaining', ...
		    n, NumExp, ListType, (t / n) * (NumExp - n) / 60);
  disp(ProgStr)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, v, I, Okay] = GetCleanData(FileName, This)
%disp(sprintf('%s, %s', FileName, This.WhichCell))
AbfS = LoadAbf(FileName);
t = AbfS.Time * 1000;  %convert to ms

FNames = fieldnames(AbfS.Units);
Current = [];
Voltage = [];
for n = 1:length(FNames)
  Unit = AbfS.Units.(FNames{n});
  if(strcmp(Unit, 'mV'))
    Voltage = [Voltage, n];
  elseif(strcmp(Unit, 'nA'))
    Current = [Current, n];
  end
end
if(length(Current) == 1 & length(Voltage == 1))
  %There are two possibilities
  %1) It's DCC, and we want the voltage that corresponds to the
  %    requested cell (top/bot).
  %2) It's two-electrode bridge mode, in which case we want the
  %    voltage trace that corresponds to the requested cell, and
  %    the current that corresponds to the other cell
  v = AbfS.Data.(FNames{Voltage});
  I = AbfS.Data.(FNames{Current});
  if(length(strfind(FNames{Voltage}, This.WhichCell)) == 0)
    Okay = false;
    return
  elseif(length(strfind(FNames{Current}, This.WhichCell)) > 0)
    Okay = true;
    DCC = true;
  else
    Okay = true;
    DCC = false;
  end
elseif(length(Current) == 2 & length(Voltage == 2))
  %There are two possibilities
  %1) It's DCC, and we want the half that corresponds to the
  %    requested cell (top/bot).
  %2) It's two-electrode bridge mode, in which case one voltage
  %    trace is good, and the other current trace corresponds to
  %    the good voltage trace; and one current trace just is a
  %    noisy zero.
  if(length(strfind(FNames{Voltage(1)}, This.WhichCell)) > 0)
    v_n = 1;
  else
    v_n = 2;
  end
  if(length(strfind(FNames{Current(1)}, This.WhichCell)) > 0)
    I_n = 1;
    IOther_n = 2;
  else
    I_n = 2;
    IOther_n = 1;
  end
  
  v = AbfS.Data.(FNames{Voltage(v_n)});
  I = AbfS.Data.(FNames{Current(I_n)});
  IOther = AbfS.Data.(FNames{Current(IOther_n)});
  if(std(IOther) < .1)
    %Bridge mode, we requested the "bad" cell
    Okay = false;
    return
  elseif(std(I) < .1)
    %Bridge mode, we requested the "good" cell
    Okay = true;
    I = IOther;
    DCC = false;
  else
    %DCC
    Okay = true;
    DCC = true;
  end
else
  error(sprintf('Problem getting data from %s', FileName))
end

if(DCC)
  [t, v, I, DCC_Info] = CleanDCC(t, v, I);
  if(~isfinite(DCC_Info.DCC_Freq))
    disp('Warning:  weird DCC signal!')
    DCC_Info
  else
    disp(sprintf('DCC Freq = %g kHz', DCC_Info.DCC_Freq))
  end
end

SmoothTime = .5;
if(t(2) - t(1) < SmoothTime)
  n = 2;
  while(t(n) - t(1) < SmoothTime)
    n = n + 2;
  end
  [B,A] = butter(2, 2 / n,'low');
  v = filtfilt(B, A, v);
end
return