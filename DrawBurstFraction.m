function DrawBurstFraction(DataStruct, ExpType, varargin)
% Usage:
%   DrawBurstFraction(Organized, 'GMGM')
%     ... OR ...
%   DrawBurstFraction(Experiments, 'ML')
% if there is no second argument, the first form is assumed.
if(nargin < 2)
  ExpType = 'GMGM';
end
if(strcmp(ExpType, 'GMGM'))
  DrawGMGMBurstFraction(DataStruct, varargin{:});
elseif(strcmp(ExpType, 'ML'))
  DrawMLBurstFraction(DataStruct, varargin{:});
else
  error(['Unknown experiment type: ', ExpType])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawGMGMBurstFraction(Organized)
Conditions = Organized.Conditions;
for n = 1:length(Conditions)
  Condition = Conditions{n};
  if(strcmp(Condition, 'ptx'))
    continue;
  end
  MakeGMGMBurstFractionFigure(Organized, Condition);
  MakeGMGMBurstFractionFigure(Organized, 'ptx', Condition);
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeGMGMBurstFractionFigure(Organized, Condition, NonControl)
g_syn = 10:15:100;
g_h = 10:15:100;

%find mean control frequency
Temp = zeros(length(g_syn), length(g_h));
Averages = Temp;
TypeArr = Organized.(sprintf('Type_%s', Condition));
if(nargin == 3)
  NonControlArray = Organized.(sprintf('Type_%s', NonControl));
  AdmitArray = isfinite(NonControlArray);
  Title = sprintf('Fraction of Bursters for %s control', NonControl);
else
  AdmitArray = ones(length(Organized.g_h), 1);
  Title = sprintf('Fraction of Bursters for %s', Condition);
end
for n = 1:length(Organized.g_h)
  gs = find(g_syn == Organized.g_syn(n));
  if(length(gs) == 0)
    continue;
  end
  gh = find(g_h == Organized.g_h(n));
  if(length(gh) == 0)
    continue;
  end
  if(~isfinite(TypeArr(n)) | ~AdmitArray(n))
    continue;
  end
  
  Temp(gs, gh) = Temp(gs, gh) + 1;
  if(TypeArr(n) > 0)
    Averages(gs, gh) = Averages(gs, gh) + 1;
  end
end
gInd = find(isfinite(Averages) & Temp > 0 & Averages > 0);
Temp = Averages ./ Temp;
[x, y] = ind2sub(size(Temp), gInd);
x = g_syn(x);
y = g_h(y);
Colors = Temp(gInd);

h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');
clf;
hold on;
CircleSize = 1000;
scatter(x, y, CircleSize, Colors, 'filled');
Map = repmat((200:-1:0)'/200, 1, 3);
colormap(Map);
colorbar;
caxis([0 1]);

xlim([0 140])
ylim([0 140])
axis square;
xlabel('g_s_y_n');
ylabel('g_h');
title(Title);
Axes_h = get(h, 'CurrentAxes');
set(Axes_h, 'FontName', 'Arial');
set(Axes_h, 'FontSize', [30]);
set(Axes_h, 'YTick', 20:20:140);
set(Axes_h, 'XTick', 0:20:140);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawMLBurstFraction(Experiments, Intrinsics)
ConditionList = GetConditionList(Experiments);
colorMapType = 'Red';

for n=1:length(ConditionList)
  CellType = ConditionList{n};
  MatchList = GetSpecifiedAnalysis(Experiments, 'CellType', CellType);
  NumCells = length(unique({MatchList.ID}));
  if(NumCells < 4)
    continue
  end
  MatchList = RemoveUnwanted(MatchList, Intrinsics);
  
  [g, Dummy, ReverseInd] = unique([cat(1, MatchList.g_syn), ...
		            cat(1, MatchList.g_h)], 'rows');
  Density = GetDensity(MatchList, ReverseInd);
  
  Title = sprintf('Burst density for %s', CellType);
  h = NamedFigure(Title);
  set(h, 'WindowStyle', 'docked');
  clf;
  hold on;
  CircleSize = 1000;
  GoodInd = find(g(:,1) < 110);
  scatter(g(GoodInd,1), g(GoodInd,2), CircleSize, Density(GoodInd), 'filled');
  switch(colorMapType)
   case 'Rainbow',
    Map = RainbowColorMap(256);
   case 'Gray',
    Map = repmat((200:-1:0)'/200, 1, 3);
   case 'Red',
    Map = repmat((200:-1:0)'/200, 1, 3);
    Map(:,1) = 1;
  end
  colormap(Map);
  colorbar;
  caxis([0 1]);
  xlim([0 115])
  ylim([0 115])
  axis square;
  xlabel('g_s_y_n', 'FontName', 'Arial', 'FontSize', 18);
  ylabel('g_h', 'FontName', 'Arial', 'FontSize', 18);
  title(Title, 'FontName', 'Arial', 'FontSize', 18);
  Axes_h = get(h, 'CurrentAxes');
  set(Axes_h, 'FontName', 'Arial');
  set(Axes_h, 'FontSize', [30]);
  set(Axes_h, 'YTick', 20:20:140);
  set(Axes_h, 'XTick', 0:20:140);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ConditionList = GetConditionList(Experiments)
ConditionList = {};
for n = 1:length(Experiments)
  Temp = Experiments(n).ConditionList;
  for m = 1:length(Temp)
    Temp{m} = StripNums(Temp{m});
  end
  ConditionList = {ConditionList{:}, Temp{:}};
end
ConditionList = unique(ConditionList);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MatchList = RemoveUnwanted(MatchList, Intrinsics)
OldMatchList = MatchList;
MatchList = [];
WantedList = {Intrinsics.ID};
for n = 1:length(OldMatchList)
  A_n = OldMatchList(n);
  if(ismember(A_n.ID, WantedList))
    MatchList = [MatchList, A_n];
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StripStr = StripNums(ConditionStr)
StripStr = regexp(ConditionStr, '[a-zA-Z]*', 'match');
if(length(StripStr) == 0)
  ErrStr = fprintf('Experiment type %s is invalid:  %s', ConditionStr, ...
		   'must contain letters');
  error(ErrStr)
else
  StripStr = StripStr{1};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Density = GetDensity(MatchList, ReverseInd)
NumPoints = max(ReverseInd);
Density = zeros(NumPoints, 1);
NumCells = 0;
for n = 1:NumPoints
  MatchInd = find(ReverseInd == n);
  NumExp = length(MatchInd);
  NumBurst = sum(cat(1, MatchList(MatchInd).Cat)==3);
  if(NumExp == 0)
    Density(n) = 0;
  else
    Density(n) = NumBurst / NumExp;
  end
end

return
