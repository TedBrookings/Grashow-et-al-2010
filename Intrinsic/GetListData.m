function varargout = GetListData(ExpList, LookIn, Subfields, ...
				 Type1, Type2)
% [List1, List2] = GetListData(ExpList, LookIn, Subfields, Type1, Type2)
%  looks through lists produeced by RunIntrinsic.m and extracts the
%  requested fields, putting results into a list.
% INPUT:
%     -ExpList:  IVList, RampList, or FIList
%     -LookIn:  A string that corresponds to the list (should be
%      'IV', 'Ramp', or 'FI' respectively)
%     -Subfields:  a cell array of strings that give the names of
%      the structure fields that contain the requested data
%     -Type1:  Experiment condition ('ptx', '5ht', 'oxo', or 'da')
%     -Type2:  OPTIONAL.  Another experimental condition.  Then
%      only experiments with both conditions will be included
% OUTPUT:
%     -List1:  A list of the values corresponding to Type1
%     -List2:  OPTIONAL.  A list of values corresponding to Type2
% EXAMPLE:
%   the following extracts spike thresholds for ptx and 5ht:
% [Thresh_ptx, Thresh_5ht] = GetListData(RampList, 'Ramp', ...
%        {'PreMaxCurve', 'V'}, 'ptx', '5ht');

TypeStr1 = ['file_', Type1];
if(nargin == 5)
  TypeStr2 = ['file_', Type2];
end

Val1 = [];
Val2 = [];
for Num = 1:length(ExpList)
  Temp1 = GetVal(ExpList(Num), LookIn, TypeStr1, Subfields);
  if(length(Temp1) == 0)
    continue;
  end
  if(nargin == 5)
    Temp2 = GetVal(ExpList(Num), LookIn, TypeStr2, Subfields);
    if(length(Temp2) == 0)
      continue;
    end
    Val2 = [Val2, Temp2];
  end
  Val1 = [Val1, Temp1];
end
if(nargin == 5)
  varargout = {Val1, Val2};
else
  varargout = {Val1};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Val = GetVal(ExpStruct, LookIn, TypeStr, Subfields)
if(length(ExpStruct.Files.(TypeStr)) == 0)
  Val = [];
  return
end
if(strcmp(class(Subfields), 'cell'))
  Val = ExpStruct.(LookIn).(TypeStr).(Subfields{1});
  for n = 2:length(Subfields)
    Val = Val.(Subfields{n});
  end
else
  Val = ExpStruct.(LookIn).(TypeStr).(Subfields);
end
return