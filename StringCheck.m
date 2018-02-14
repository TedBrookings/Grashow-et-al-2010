function IsMatch = StringCheck(CheckString, Pattern, StrictCase)
% IsMatch = StringCheck(CheckString, Pattern, StrictCase)
% Checks to see if CheckString contains Pattern.
%  INPUT PARAMETERS:
%   -CheckString: string to be checked
%   -Pattern: pattern to test against
%    OPTIONAL:
%    -StrictCase:  If set to true, upper/lower case must be exact.
%           (Defaults to false)
%  OUTPUT PARAMETERS:
%   -IsMatch:  boolean specifying if there is a match
if(nargin < 3)
  StrictCase = false;
end

if(~StrictCase)
  CheckString = lower(CheckString);
  Pattern = lower(Pattern);
end

IsMatch = (length(strfind(CheckString, Pattern)) > 0);
return