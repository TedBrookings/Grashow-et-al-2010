function DrawScatter(Organized)

Conditions = Organized.Conditions;
for n = 1:length(Conditions)
  Condition = Conditions{n};
  if(strcmp(Condition, 'ptx'))
    continue;
  end

  X = Organized.Type_ptx;
  TypeString = sprintf('Type_%s', Condition);
  Y = Organized.(TypeString);
  MakeFreqHist(Y, Condition);
  MakeFreqHist(X, Condition, Y);
  
  IndAltSpikes = MakeShiftHist(X, Y, Condition, Organized);
  
  MakeScatterFig(X, Y, Condition, IndAltSpikes);
  MakeCloseUpFig(X, Y, Condition, IndAltSpikes);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeScatterFig(X, Y, Condition, IndAltSpikes)
Title = sprintf('ptx vs %s', Condition);
h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');

NumDiv = ceil(max(X) / .1);
Colors = GetColors(NumDiv);

XStop = NumDiv * .1;
YStop = .1 * round(max(Y) / .1);

clf;
plot([0,XStop], [0, XStop], 'k--');
hold on;
Ind = find((X == 0 | Y == 0) & isfinite(X) & isfinite(Y));
Ind1 = setdiff(Ind, IndAltSpikes);
Ind2 = intersect(Ind, IndAltSpikes);
plot(X(Ind1), Y(Ind1), 'bo', 'MarkerSize', 7);
plot(X(Ind2), Y(Ind2), 'b*', 'MarkerSize', 10);
NumDots = length(Ind);
for Div = 1:NumDiv
  if(Div == 0)
    Ind = find(0 < X & X <= .1);
  else
    Ind = find((Div-1) * .1 < X & X <= Div * .1 & Y > 0);
  end
  Ind1 = setdiff(Ind, IndAltSpikes);
  Ind2 = intersect(Ind, IndAltSpikes);
  
  NumDots = NumDots + length(Ind);
  if(Div ~= NumDiv)
    plot([Div * .1, Div * .1], [0, YStop], ':', 'Color', [.5,.5,.5]);
  end
  if(length(Ind) == 0)
      continue;
  end
  plot(X(Ind1), Y(Ind1), 'o', 'MarkerEdgeColor', 'k', ...
       'MarkerSize', 7, 'MarkerFaceColor', Colors(Div,:));
  plot(X(Ind2), Y(Ind2), '*', 'MarkerEdgeColor', Colors(Div,:), ...
       'MarkerSize', 10);
end
hold off;

title(Title);
xlabel('ptx Frequency (Hz)');
ylabel(sprintf('%s Frequency (Hz)', Condition));
xlim([0 XStop]);
ylim([0 YStop]);
Axes_h = get(h, 'CurrentAxes');
set(Axes_h, 'FontName', 'Ariel');
set(Axes_h, 'FontSize', [15]);
set(Axes_h, 'XTick', 0:.1:XStop);
disp(sprintf('Number of dots for scatter plot of %s:  %g', ...
	     Title, NumDots))
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeCloseUpFig(X, Y, Condition, IndAltSpikes)
Title = sprintf('ptx vs %s close-up', Condition);
h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');

NumDiv = ceil(max(X) / .1);
Colors = GetColors(NumDiv);

XStop = .3;
YStop = .4;

clf;
plot([0, XStop], [0, XStop], 'k--');
hold on;
Div = 2;
Ind = find((Div-1) * .1 < X & X <= Div * .1 & Y ~= 0);
Ind1 = setdiff(Ind, IndAltSpikes);
Ind2 = intersect(Ind, IndAltSpikes);
plot(X(Ind1), Y(Ind1), 'o', 'MarkerEdgeColor', 'k', ...
     'MarkerSize', 14,'MarkerFaceColor', Colors(Div,:));
plot(X(Ind2), Y(Ind2), '*', 'MarkerEdgeColor', Colors(Div,:), ...
     'MarkerSize', 20);
hold off;

title(Title);
xlabel('ptx Frequency (Hz)');
ylabel(sprintf('%s Frequency (Hz)', Condition));
xlim([0 .3]);
ylim([0 .4]);
AxesHandle = get(h, 'CurrentAxes');
set(AxesHandle, 'XTick', 0:.05:XStop);
set(AxesHandle, 'YTick', 0:.05:YStop);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Colors = GetColors(NumDiv)
Colors = [[1 0 0]; ...
	  [0 1 0]; ...
	  [0 0 1]; ...
	  [1 .5 0]; ...
	  [0 .5 1]; ...
	  [.5 1 0]; ...
	  [.5 0 1]; ...
	  [0 1 .5]; ...
	  [1 0 .5]; ...
	  [0 1 1]; ...
	  [1 1 0]; ...
	  [1 0 1]; ...
	  [0 0 0]; ...
	  [1 1 1]];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeFreqHist(Freq, Condition, varargin)
if(nargin == 2)
  Title = sprintf('Histogram of %s frequencies', Condition);
else
  Freq2 = varargin{1};
  GoodInd = find(isfinite(Freq) & isfinite(Freq2));
  Freq = Freq(GoodInd);
  Title = sprintf('Histogram of ptx (vs. %s)', Condition);
end
h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');
Bins = [0:.1:(max(Freq) + .1)];
hist(Freq, Bins);
xlabel('Freq (Hz)');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IndAltSpikes = MakeShiftHist(X, Y, Condition, Organized)

SpikesPerBurst_X = GetSpikesPerBurst(Organized.Analysis_ptx);
AnalysisString = sprintf('Analysis_%s', Condition);
SpikesPerBurst_Y = GetSpikesPerBurst(Organized.(AnalysisString));
IndAltSpikes = find(SpikesPerBurst_X == 1 | SpikesPerBurst_Y == 1);
IndBurst = find(SpikesPerBurst_X > 1 & SpikesPerBurst_Y > 1);
[Sig, PVal] = ttest(X, Y, .05, 'both');
MessageOut = sprintf('T-Test for frequency, ptx vs %s: p = %g', ...
		     Condition, PVal);
disp(MessageOut)
[Sig, PVal] = ttest(X(IndAltSpikes), Y(IndAltSpikes), .05, 'both');
MessageOut = sprintf('T-Test for frequency, ptx vs %s (%s): p = %g', ...
		     Condition, 'one spike per burst', PVal);
disp(MessageOut)
[Sig, PVal] = ttest(X(IndBurst), Y(IndBurst), .05, 'both');
MessageOut = sprintf('T-Test for frequency, ptx vs %s (%s): p = %g', ...
		     Condition, ' >1 spike per burst', PVal);
disp(MessageOut)

Title = sprintf('Histogram of shifts: ptx vs %s', Condition);
h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');

Ind = find(X > 0 & Y > 0);
Shifts = (Y(Ind) - X(Ind));
[N, XBin] = hist(Shifts, 100);
hist(Shifts, 100);
[NMax, MaxInd] = max(N);
XMax = XBin(MaxInd);
FirstAbove = find(XBin > XMax & N < exp(-9/2), 1);
XAbove = XBin(FirstAbove);
MessageOut = sprintf('In %s, Large Shift = %g', Condition, XAbove);
disp(MessageOut)
hold on
plot([XAbove, XAbove], [0, NMax], 'r-')
xlabel(sprintf('%s Freq - ptx Freq (Hz)', Condition))
hold off

Ind2 = find(Shifts >= XAbove);
Ind2 = Ind(Ind2);
Ind1 = find(Shifts < XAbove);
Ind1 = Ind(Ind1);
Title = sprintf('Spikes Per Burst for Large Shifts (%s)', Condition);
h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');
plot(SpikesPerBurst_X(Ind1), SpikesPerBurst_Y(Ind1), ...
     'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b')
hold on
plot(SpikesPerBurst_X(Ind2), SpikesPerBurst_Y(Ind2), ...
     'r*', 'MarkerSize', 12, 'LineWidth', 1.5)
hold off
xlabel('Mean spikes per burst in ptx', 'FontSize', 16);
ylabel(sprintf('Mean spikes per burst in %s', Condition), 'FontSize', 16);

Title = sprintf('Spikes Per Burst Hist for Large Shifts (ptx vs %s)', Condition);
h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');
hist(SpikesPerBurst_X(Ind2), 10);
Title = sprintf('Spikes Per Burst Hist for Large Shifts (%s)', Condition);
h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');
hist(SpikesPerBurst_Y(Ind2), 10);

return
for n=1:length(Ind)
  PrintWeird(Organized, Ind(n), Condition);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SpikesPerBurst = GetSpikesPerBurst(AnalysisList)
Len = length(AnalysisList);
SpikesPerBurst = zeros(Len, 2);
for n=1:Len
  Temp = AnalysisList(n);
  if(length(Temp.Cell0) == 0 | length(Temp.Cell1) == 0 ...
     | Temp.HalfCenter.Freq == 0)
    continue
  end
  try
    SpikesPerBurst(n,1) = Temp.Cell0.Burst.NumSpikes.Mean;
    SpikesPerBurst(n,2) = Temp.Cell1.Burst.NumSpikes.Mean;
    if(SpikesPerBurst(n,1) == 0)
      SpikesPerBurst(n,1) = 1;
    end
    if(SpikesPerBurst(n,2) == 0)
      SpikesPerBurst(n,2) = 1;
    end
  catch
    Temp
    whos
    error('poo')
  end
end

SpikesPerBurst = max(SpikesPerBurst, [], 2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintWeird(Organized, n, Condition)
if(ispc)
  Slash = '\';
else
  Slash = '/';
end

g_syn = Organized.Analysis_ptx(n).g_syn;
g_h = Organized.Analysis_ptx(n).g_h;
ExpNum = Organized.Analysis_ptx(n).ExpNum;
ptxFile = Organized.Analysis_ptx(n).FileName;
CondFile = Organized.(sprintf('Analysis_%s', Condition))(n).FileName;

Ind = strfind(ptxFile, Slash);
if(length(Ind) > 0)
  ptxFile = ptxFile((Ind(end)+1):end);
end
Ind = strfind(CondFile, Slash);
if(length(Ind) > 0)
  CondFile = CondFile((Ind(end)+1):end);
end

OutLine = sprintf('%s %s Exp# %g g_syn %g g_h %g', ptxFile, CondFile, ...
		  ExpNum, g_syn, g_h);
disp(OutLine)
return