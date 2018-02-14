function options = GetOptions(defaults, args, ignoreExtra)
%options = GetOptions(defaults, args)
% INPUTS:
%  -defaults: a cell array with keyword, value ordering
%             or a structure (although that's rather inconvenient)
%  OPTIONAL
%  -args: either a structure or
%                a cell array
%         if a cell array, user can pass ALL keywords or
%                                   just go by input order.
%         if args is ommitted, the options will be the defaults.
%  -ignoreExtra: if set to true, ignore options specified in args
%                  that don't correspond to anything in defaults.
  %                This defaults to false
% OUTPUTS:
%  -options: structure with fields set to their proper values

if nargin < 1 || nargin > 3
  help GetOptions
  error('Invalid number of inputs.')
end
if ~iscell(defaults) && ~isstruct(defaults)
  help GetOptions
  error('"defaults" must be a cell array or structure.')
end
if iscell(defaults);
  defaults = cellOptionsToStruct(defaults);
end

if nargin == 1 || length(args) == 0
  options = defaults;
  return
end

if nargin < 3
  ignoreExtra = false;
end

if isstruct(args)
  options = getStructOptions(defaults, args, ignoreExtra);
elseif iscell(args)
  options = getCellArrOptions(defaults, args, ignoreExtra);
else
  help GetOptions
  error('"args" must be a structure or a cell array.')
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optStruct = cellOptionsToStruct(cellOpts)
numOpts = length(cellOpts);
if mod(numOpts, 2) ~= 0
  error('options must come in "key" "value" pairs.')
end
numOpts = numOpts / 2;
for n = 1:numOpts
  mVal = 2 * n;
  mKey = mVal - 1;
  key = cellOpts{mKey};
  val = cellOpts{mVal};
  if ~ischar(key)
    disp(key)
    error('Invalid option key.\n')
  end
  optStruct.(key) = val;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function options = getStructOptions(defaults, args, ignoreExtra)
options = defaults;
fNames = fieldnames(args);
for n = 1:length(fNames)
  field_n = fNames{n};
  if isfield(options, field_n)
    options.(field_n) = args.(field_n);
  elseif ~ignoreExtra
    error('%s is an invalid option', field_n)
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function options = getCellArrOptions(defaults, args, ignoreExtra)
fNames = fieldnames(defaults);
n = 1;
pos = 1;
argStruct = [];
while n < length(args)
  if isfield(defaults, args{n})
    argStruct.(args{n}) = args{n+1};
    n = n + 2;
  else
    if pos <= length(fNames)
      argStruct.(fNames{pos}) = args{n};
      pos = pos + 1;
      n = n + 1;      
    elseif ~ignoreExtra
      error('Too many arguments.\n')
    end
  end
end
options = getStructOptions(defaults, argStruct, ignoreExtra);
return