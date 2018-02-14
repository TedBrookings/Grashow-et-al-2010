function DrawMap(Map, varargin)
%varargin is ValueString in GetMapValue()
if(ischar(Map))
  load(Map);   %if Map is a string, assume it's the filename to the
               %map file
end

if(nargin == 2)
  Title = sprintf('3D map for %s', Map.Type);
else
  Title = sprintf('3D map for %s (with ptx)', Map.Type);
end

h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');
clf;

%Check to see if this is a integer map
IsInt = IsIntegerMap(Map);
if(IsInt & length(Map.Ind) > 0)
  Values = Map.Value(Map.Ind);
  scatter(Map.g_syn(Map.Ind), Map.g_h(Map.Ind), 100, Values, 'filled');
  CustomMap = MakeCustomColorMap(Map.Values);
  colormap(CustomMap);
  if(max(Values > min(Values)))
    colorbar;
  end
else
  if(nargin == 3)
    Make3DFreqPlot(Map, varargin{2});
  else
    Make3DFreqPlot(Map);
  end
end

if(sum(Map.g_syn(Map.Ind) > 120) + sum(Map.g_h(Map.Ind) > 120) == 0)
  axis([0 120 0 120]);
end
xlabel('g_s_y_n (nS)','FontName', 'Arial', 'FontSize', [30]);
ylabel('g_h (nS)','FontName', 'Arial', 'FontSize', [30]);
zlabel('Frequency (Hz)','FontName', 'Arial', 'FontSize', [30]);
title(Title);
Axes_h = get(h, 'CurrentAxes');
set(Axes_h, 'FontName', 'Arial');
set(Axes_h, 'FontSize', [20]);
set(Axes_h, 'GridLineStyle', ':','LineWidth', 2);

if(~IsInt)
  view(3);
  grid on;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IsInt = IsIntegerMap(Map)
Values = Map.Values(Map.Ind);
if(sum(round(Values) ~= Values) == 0)
  IsInt = 1;
else
  IsInt = 0;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CustomMap = MakeCustomColorMap(Values)
Map = [[0 0 0]; ...
       [0 1 0]; ...
       [0 0 1];
       [1 0 0]];
MaxColors = size(Map, 1) - 1;

%Check if each value occurs
CustomMap = [];
for n = 0:MaxColors
  if(length(find(Values == n)) > 0)
    CustomMap = [CustomMap; Map(n + 1,:)];
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Make3DFreqPlot(Map, Map2)
NumPoints = length(Map.Cats);
BadInd = setdiff(1:NumPoints, Map.Ind);

%plot3(x, y, z, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 15);
hold on;
for n = 1:length(Map.Ind)
  m = Map.Ind(n);
  x = Map.g_syn(m);
  y = Map.g_h(m);
  z = Map.Values(m);
  switch(Map.Cats(m))
   case 0, C_m = 'k';
   case 1, C_m = 'g';
   case 2, C_m = 'b';
   case 3, C_m = 'r';
   otherwise, error('Invalid Category');
  end
  plot3(x, y, z, [C_m, 'o'], 'MarkerFaceColor', C_m, 'MarkerSize', 16);
  plot3([x, x], [y, y], [0, z], [C_m, '-'],'LineWidth', 2);
  if(nargin == 2)
    k = find(Map2.g_syn == x & Map2.g_h == y & Map2.Cats == 3);
    if(length(k) > 0)
      z2 =  Map2.Values(k);
      plot3([x, x], [y, y], [0, z2], 'k--','LineWidth', 2);
      plot3(x, y, z2, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 13);
      plot3(x, y, z2, 'kd', 'MarkerFaceColor', 'r', 'MarkerSize', 10); 
    end
  end
end

if(nargin == 2)
  for n = 1:length(Map2.Ind)
    m = Map2.Ind(n);
    x = Map2.g_syn(m);
    y = Map2.g_h(m);
    z = Map2.Values(m);
    switch(Map2.Cats(m))
     case 0, C_m = 'k';
     case 1, C_m = 'g';
     case 2, C_m = 'b';
     case 3, C_m = 'r';
     otherwise, error('Invalid Category');
    end
    plot3([x, x], [y, y], [0, z], 'k--','LineWidth', 2);
    
    plot3(x, y, z, 'kd', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    plot3(x, y, z, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 13);
  end
else
  for n = 1:length(BadInd)
    m = BadInd(n);
    x = Map.g_syn(m);
    y = Map.g_h(m);
    z = 0;
    switch(Map.Cats(m))
     case 0, C_m = 'k';
     case 1, C_m = 'g';
     case 2, C_m = 'b';
     case 3, C_m = 'r';
     case 4, C_m = 'r';   %remove this, it shouldn't be there!
     otherwise, error('Invalid Category');
    end
    plot3(x, y, z, [C_m, 'o'], 'MarkerFaceColor', C_m, 'MarkerSize', 13);
  end
end
hold off;
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
