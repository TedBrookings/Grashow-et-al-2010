function DrawDecrease(ControlMap, ModMap, DirName, PlotTrendBucking)
Ind = sort([strfind(DirName, '\'), strfind(DirName, '/')]);
TitleBase = DirName((Ind(end)+1):end);
Title = sprintf('%s: %s', TitleBase, ModMap.Type);
Ind = strfind(Title, '_');
Title(Ind) = '-';

MinusChangeThreshold = .001;  %absolute threshold for negative change
switch(ModMap.Type)
 case '5ht',
  LargeShift = .345;
 case 'oxo',
  LargeShift = .125;
 case 'DC',
  LargeShift = .2;
 otherwise,
  error(sprintf('LargeShift not defined for %s', ModMap.Type))
end
PlusChangeThreshold = LargeShift;    %threshold for z-scored positive change

ArrowScale = 1.5;
CircleSize = 800;
if(nargin < 4)
  PlotTrendBucking = '';
end
if(strcmp(PlotTrendBucking, 'ptx'))
  Bucks = GetBucks(ControlMap);
  BucksLineStyle = 'g-';
  h = NamedFigure([Title, ' with ptx bucks']);
  set(h, 'WindowStyle', 'docked');
  clf;
  hold on;
  Ind = find(ControlMap.Cats == 3);
  scatter(ControlMap.g_syn(Ind), ControlMap.g_h(Ind), ...
	  CircleSize, 'r', 'filled');
elseif(strcmp(PlotTrendBucking, 'mod'))
  Bucks = GetBucks(ModMap);
  BucksLineStyle = 'b-';
  h = NamedFigure(sprintf('%s with %s bucks', Title, ModMap.Type));
  set(h, 'WindowStyle', 'docked');
  clf;
  hold on;
  Ind = find(ModMap.Cats == 3);
  scatter(ModMap.g_syn(Ind), ModMap.g_h(Ind), ...
	  CircleSize, 'r', 'filled');
else
  h = NamedFigure(Title);
  set(h, 'WindowStyle', 'docked');
  clf;
  hold on;
end


%ControlMap.g_syn = ControlMap.g_syn(ControlMap.Ind);
%ControlMap.g_h = ControlMap.g_h(ControlMap.Ind);
%ControlMap.Values = ControlMap.Values(ControlMap.Ind);

%ModMap.g_syn = ModMap.g_syn(ModMap.Ind);
%ModMap.g_h = ModMap.g_h(ModMap.Ind);
%ModMap.Values = ModMap.Values(ModMap.Ind);

%Find the list of networks that exist in both conditions
[g, IndC, IndM] = intersect([ControlMap.g_syn' ControlMap.g_h'], ...
			    [ModMap.g_syn' ModMap.g_h'], 'rows');
if(length(IndC) == 0)
  disp(sprintf('Warning, no overlapping networks for %s', ModMap.Type))
  return
end

ModCat = ModMap.Cats(IndM);
ModVal = ModMap.Values(IndM);
ModVal(find(ModCat ~= 3)) = 0;
ControlCat = ControlMap.Cats(IndC);
ControlVal = ControlMap.Values(IndC);
ControlVal(find(ControlCat ~= 3)) = 0;
Change = ModVal - ControlVal;
zChange = Change ./ (ModVal + ControlVal);
zChange = zscore(zChange);

for n = 1:length(Change)
  if(ModVal(n) == 0 | ControlVal(n) == 0)
    continue;
  end
  DrawChange(g(n,1), g(n,2), Change(n), zChange(n), MinusChangeThreshold, ...
	     PlusChangeThreshold, ArrowScale, CircleSize);
end

if(length(PlotTrendBucking) > 0)
  DrawBucks(Bucks, BucksLineStyle);
end

hold off;
title(Title,'FontName', 'Arial', 'FontSize', [40]);
xlim([0 140]);
ylim([0 140]);
xlabel('g_s_y_n', 'FontName', 'Arial', 'FontSize', [40]);
ylabel('g_h','FontName', 'Arial', 'FontSize', [40]);
Axes_h = get(h, 'CurrentAxes');
set(Axes_h, 'FontName', 'Arial');
set(Axes_h, 'FontSize', [30]);
axis square
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawChange(g_syn, g_h, Change, zChange, MinusChangeThreshold, ...
		    PlusChangeThreshold, ArrowScale, CircleSize)
scatter(g_syn, g_h, CircleSize, 'r', 'filled');

Scale = ArrowScale;
if(Change > PlusChangeThreshold)
  %XPoints = [g_syn - 2 * Scale, g_syn, g_syn + 2 * Scale];
  %YPoints = [g_h + 3 * Scale, g_h + 5 * Scale, g_h + 3 * Scale];
  %Code = 'g-';
  XPoints = g_syn;
  YPoints = g_h;
  Width = 3;
  Code = {'k*', 'LineWidth', Width, 'MarkerSize', 40};
elseif(Change < -MinusChangeThreshold)
  XPoints = [g_syn - 2 * Scale, g_syn, g_syn + 2 * Scale];
  YPoints = [g_h - 3 * Scale, g_h - 5 * Scale, g_h - 3 * Scale];
  
  Width = 6;
  Code = {'k-', 'LineWidth', Width};
  plot([g_syn, g_syn], [g_h - 5 * Scale, g_h + 5 * Scale], Code{:});
else
  return
end
plot(XPoints, YPoints, Code{:});
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bucks = GetBucks(ControlMap)
g_h = ControlMap.g_h;
g_syn = ControlMap.g_syn;
vals = ControlMap.Values;
ZeroInd = find(ControlMap.Cats ~= 3);
vals(ZeroInd) = 0;

delta_h = GetIncrement(g_h);
delta_syn = GetIncrement(g_syn);

Bucks = [];

NumExp = length(g_h);
for n = 1:(NumExp-1)
  for m = (n+1):NumExp
    d_h = g_h(n) - g_h(m);
    d_syn = g_syn(n) - g_syn(m);
    
    if(abs(d_h) == delta_h & d_syn == 0)
      grad_h = (vals(n) - vals(m)) / d_h;
      if(grad_h < 0)
	Buck.Type = 'h';
	Buck.g_syn_1 = g_syn(n);
	Buck.g_syn_2 = g_syn(m);
	Buck.g_h_1 = g_h(n);
	Buck.g_h_2 = g_h(m);
	if(IsExpected(g_syn, g_h, vals, n, m, Buck))
	    continue;
	end
	Bucks = [Bucks, Buck];
      end
    elseif(abs(d_syn) == delta_syn & d_h == 0)
      grad_syn = (vals(n) - vals(m)) / d_syn;
      if(grad_syn > 0)
	Buck.Type = 'syn';
	Buck.g_syn_1 = g_syn(n);
	Buck.g_syn_2 = g_syn(m);
	Buck.g_h_1 = g_h(n);
	Buck.g_h_2 = g_h(m);
	if(IsExpected(g_syn, g_h, vals, n, m, Buck))
	    continue;
	end
	Bucks = [Bucks, Buck];
      end
    end
   end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Expected  = IsExpected(g_syn, g_h, vals, n, m, Buck)
Expected = false;
if(vals(n) > 0)
  if(vals(m) > 0)
    return
  end
else
  Temp = m;
  m = n;
  n = m;
end 
%now we know that vals(m) = 0, vals(n) > 0
if(strcmp(Buck.Type, 'h'))
  %know g_h(m) > g_h(n), otherwise wouldn't be a trend-buck
  Ind = find(g_h > g_h(m) & g_syn == g_syn(m) & vals > 0);
  if(length(Ind) == 0)
    Expected = true;
  end
else
  %know g_syn(m) < g_syn(n), otherwise wouldn't be a trend-buck
  Ind = find(g_syn < g_syn(m) & g_h == g_h(m) & vals > 0);
  if(length(Ind) == 0)
    Expected = true;
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta = GetIncrement(Arr)
Ind = find(Arr ~= Arr(1));
delta = min(Arr(Ind) - Arr(1));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawBucks(Bucks, LineStyle)
Width = 6;
for n=1:length(Bucks)
  Buck = Bucks(n);
  plot([Buck.g_syn_1, Buck.g_syn_2], [Buck.g_h_1, Buck.g_h_2], ...
       LineStyle, 'LineWidth', Width);
end
return
