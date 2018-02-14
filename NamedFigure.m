function varargout = NamedFigure(Name)
% h = NamedFigure(Name)
%  Looks for a figure with name Name, sets it to the current
%  figure, and returns a handle h to the figure.  If no such figure
%  exists, it creates one, names it, and returns a handle
%    INPUT:
%     -Name   The name of the figure
%    OUTPUT:
%     -h      Handle to the figure

h = findobj('Name', Name);
if(length(h) == 0)
  set(0, 'defaultaxesfontname', 'Arial')
  set(0, 'defaulttextfontname', 'Arial')
  set(0, 'defaultaxesfontsize', 15)
  set(0, 'defaulttextfontsize', 18)  
  
  h = figure;
  set(h, 'Name', Name)
  set(h, 'Renderer', 'painters')
  set(h, 'RendererMode', 'manual')
%  title(Name);  %No point to this, since it gets overwritten.
elseif(ishandle(h))
  set(0, 'CurrentFigure', h)
else
  h = figure(h);
end

switch(nargout)
 case 0, varargout = {};
 case 1, varargout = {h};
 otherwise, error('Too many output arguments.');
end
return