function varargout = RainbowColorMap(NumColors)
% CMap = RainbowColorMap(NumColors)
%
% Creates a colormap, similar to jet, but starting with black (as
% opposed to dark blue, and ending with bright red (as opposed to
% dark red).  If the colormap is not returned, it is set to be the
% active colormap.
%  INPUT PARAMETERS: (OPTIONAL)
%   -NumColors:  number of colors in colormap.  Defaults to 256.
%  OUTPUT PARAMETERS: (OPTIONAL)
%   -CMap:  the colormap matrix.  If CMap is not returned,
%     colormap(CMap) is called to set the active colormap.

if(nargin < 1)
  NumColors = 256;
end

CMap = zeros(NumColors, 3);

NumPhases = 3;

Len = round(.5 * NumColors / NumPhases);
Start_n = 2;
Stop_n = Len;
for n = Start_n:Stop_n
  CMap(n,3) = (n - 1) / (Len - 1);
end
NumPhases = NumPhases - .5;
NumColors = NumColors - Len;

Len = round(NumColors / NumPhases);
Start_n = Stop_n + 1;
Stop_n = Start_n + Len - 1;
for n = Start_n:Stop_n
  CMap(n,3) = 1;
  CMap(n,2) = (n - Start_n) / (Len - 1);
end
NumPhases = NumPhases - 1;
NumColors = NumColors - Len;

Len = round(NumColors / NumPhases);
Start_n = Stop_n + 1;
Stop_n = Start_n + Len - 1;
for n = Start_n:Stop_n
  CMap(n,3) = (Stop_n - n) / (Len - 1);
  CMap(n,2) = 1;
  CMap(n,1) = (n - Start_n) / (Len - 1);
end
NumPhases = NumPhases - 1;
NumColors = NumColors - Len;

Len = NumColors;
Start_n = Stop_n + 1;
Stop_n = Start_n + Len - 1;
for n = Start_n:Stop_n
  CMap(n,2) = (Stop_n - n) / (Len - 1);
  CMap(n,1) = 1;
end

if(nargout >= 1)
  varargout = {CMap};
else
  varargout = {};
  colormap(CMap);
end
return