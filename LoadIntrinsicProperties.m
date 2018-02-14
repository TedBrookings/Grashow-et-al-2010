function Props = LoadIntrinsicProperties
if(ispc)
  FileName = ['C:\Documents and Settings\Rachel\My Documents\', ...
	      'matlab-data-analysis\dynamic_clamp\IP summary for ted.txt'];
else
  FileName = '/home/ted/data/rachel/IP summary for ted.txt';
end
fid = fopen(FileName);
if(fid < 0)
  error(sprintf('File %s does not exist', FileName))
end

Props = [];

TextLine = fgets(fid);  %header line
TextLine = fgets(fid);
while(ischar(TextLine))
  try
    PropRec = GetPropRecord(TextLine);
  catch
    fclose(fid);
    error(sprintf('Error reading %s', FileName));
  end
  Props = [Props, PropRec];
  TextLine = fgets(fid);
end
fclose(fid);

[ExpNum, Ind] = sort(cat(1,Props.ExpNum));
Props = Props(Ind);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PropRec = GetPropRecord(TextLine)
Columns = {{'Book', '%g'}, {'ExpNum', '%g'}, {'WhichCell', '%s'}, ...
	   {'Modulator', '%s'}, {'R_ptx', '%g'}, {'R_mod', '%g'}, ...
	   {'VRest_ptx', '%g'}, {'VRest_mod', '%g'}, ...
	   {'VSpike_ptx', '%g'}, {'VSpike_mod', '%g'}};
PropRec = [];
for n = 1:length(Columns)
  [PropRec, TextLine] = GetColInfo(PropRec, TextLine, Columns{n});
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PropRec, TextLine] = GetColInfo(PropRec, TextLine, Column)
[PropRec.(Column{1}), Count, ErrMsg, NextInd] ...
    = sscanf(TextLine, Column{2}, 1);
TextLine = TextLine(NextInd:end);
return