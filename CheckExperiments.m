function CheckExperiments(ListFile, BaseDir)
% CheckExperiments(ListFile, BaseDir)
% Checks 'Experiment' files for errors
%    -BaseDir is the directory for all the experiment foloders
%         (defaults to C:\Documents and Settings\Rachel\My Documents\spectra)
%    -ListFile is a text file listing all the directories to open
%         (defaults to folder_names.txt)

if(nargin > 3)
  error('CheckExperiments.m requires 2 or fewer input arguments');
end
if(nargin < 2)
  if(ispc)
    BaseDir = 'C:\Documents and Settings\Rachel\My Documents\spectra\';
    AltDir = '/mnt/dwidget/';
  else
    BaseDir = '/mnt/dwidget/';
    AltDir = 'C:\Documents and Settings\Rachel\My Documents\spectra\';
  end
end
if(nargin < 1)
  ListFile = [BaseDir, 'folder_names.txt'];
end


[DirList, ExpList] = GetDirList(BaseDir, ListFile);
CheckExperimentFiles(DirList, ExpList, BaseDir, AltDir);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DirList, ExpList] = GetDirList(BaseDir, ListFile)
DirList = {};
ExpList = {};
fid = fopen(ListFile, 'r');
NextLine = fgetl(fid);
while ischar(NextLine);
  ExpList = {ExpList{:}, NextLine};
  DirList = {DirList{:}, [BaseDir, NextLine]};
  NextLine = fgetl(fid);
end
fclose(fid);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CheckExperimentFiles(DirList, ExpList, BaseDir, AltDir)
if(ispc)
  Slash = '\';
  AltSlash = '/';
else
  Slash = '/';
  AltSlash = '\';
end

for n = 1:length(DirList)
  Exp = ExpList{n};
  AltFileName = sprintf('%s%sExperiment_%s.txt', DirList{n}, Slash, Exp);
  Ind = strfind(Exp, '_');
  Exp = Exp((Ind+1):end);
  FileName = sprintf('%s%sExperiment_%s.txt', DirList{n}, Slash, Exp);
  CheckFile(FileName, AltFileName, ExpList{n});
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CheckFile(FileName, AltFileName, ExpNum)
gsyn = [];
gh = [];
FileNum = [];
Cond = {};
fid = fopen(FileName, 'r');
if(fid < 0)
  fid = fopen(AltFileName, 'r');
  if(fid > 0)
    ErrStr = sprintf('Rename %s to %s', StripFileName(AltFileName), ...
		     StripFileName(FileName));
    %Technically, this isn't an "error":
    %disp(ErrStr)
  else
    ErrStr = sprintf('Couldn''t open %s or %s', StripFileName(AltFileName), ...
		     StripFileName(FileName));
    disp(ErrStr)
    return
  end
end
NextLine = fgetl(fid);  %first line is a dummy
NextLine = fgetl(fid);

n = 1;
while ischar(NextLine);
  [gsyn(n), Count, ERRMSG, NextInd] = sscanf(NextLine,'%g', 1);
  NextLine = NextLine(NextInd:end);
  [gh(n), Count, ERRMSG, NextInd] = sscanf(NextLine,'%g', 1);
  NextLine = NextLine(NextInd:end);
  [FileNum(n), Count, ERRMSG, NextInd] = sscanf(NextLine,'%g', 1);
  NextLine = NextLine(NextInd:end);
  [Cond{n}, Count, ERRMSG, NextInd] = sscanf(NextLine,'%s', 1);
  NextLine = NextLine(NextInd:end);
  
  NextLine = fgetl(fid);
  n = n + 1;
end
fclose(fid);

Ind = find(strcmp(Cond, 'ptx'));
gsyn_ptx = gsyn(Ind);
gh_ptx = gh(Ind);
FileNum_ptx = FileNum(Ind);

CondList = unique(Cond);
for m = 1:length(CondList)
  Ind = find(strcmp(Cond, CondList{m}));
  Cond_m = CondList{m};
  gsyn_m = gsyn(Ind);
  gh_m = gh(Ind);
  FileNum_m = FileNum(Ind);
  CheckCondition(ExpNum, Cond_m, gsyn_m, gh_m, FileNum_m);
  if(~strcmp(Cond_m, 'ptx'))
    CheckAcross(ExpNum, ...
		'ptx', gsyn_ptx, gh_ptx, FileNum_ptx, ...
		Cond_m, gsyn_m, gh_m, FileNum_m);
    CheckAcross(ExpNum, ...
		Cond_m, gsyn_m, gh_m, FileNum_m, ...
    		'ptx', gsyn_ptx, gh_ptx, FileNum_ptx);
  end
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StripName = StripFileName(FileName)
if(ispc)
  Slash = '\';
  AltSlash = '/';
else
  Slash = '/';
  AltSlash = '\';
end
Ind = strfind(FileName, Slash);
if(length(Ind) > 0)
  StripName = FileName((Ind(end)+1):end);
else
  StripName = FileName;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CheckCondition(ExpNum, Condition, gsyn, gh, FileNum)
ErrPre = sprintf('In %s %s:', ExpNum, Condition);
[FileNum, Ind] = sort(FileNum);
gsyn = gsyn(Ind);
gh = gh(Ind);
NumLines = length(FileNum);
LastNum = FileNum(1);
LastSyn = gsyn(1);
LastH = gh(1);
for n=2:NumLines
  if(LastNum == FileNum(n))
    ErrLine = sprintf('%s FileNum=%g used by gsyn/gh = %g/%g and %g/%g', ...
		      ErrPre, LastNum, LastSyn, LastH, gsyn(n), gh(n));
    disp(ErrLine)
  end
  LastNum = FileNum(n);
  LastSyn = gsyn(n);
  LastH = gh(n);
end

UniqueSyn = unique(gsyn);
UniqueH = unique(gh);
for n=1:length(UniqueSyn)
  for m=1:length(UniqueH)
    Ind = find(gsyn == UniqueSyn(n) & gh == UniqueH(m));
    if(length(Ind) == 0)
         ErrLine = sprintf('%s missing gsyn/gh = %g/%g', ...
			   ErrPre, UniqueSyn(n), UniqueH(m));
	 disp(ErrLine)
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CheckAcross(ExpNum, ...
		     Cond_ptx, gsyn_ptx, gh_ptx, FileNum_ptx, ...
		     Cond_m, gsyn_m, gh_m, FileNum_m)

for n = 1:length(gsyn_ptx)
  Ind = find(gsyn_m == gsyn_ptx(n) & gh_m == gh_ptx(n));
  if(length(Ind) == 0)
    ErrStr = sprintf('In %s: gsyn/gh = %g/%g present in %s but not %s', ...
		     ExpNum, gsyn_ptx(n), gh_ptx(n), Cond_ptx, Cond_m);
    disp(ErrStr)
  end
end

return