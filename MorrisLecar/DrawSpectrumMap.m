function DrawSpectrumMap(Map, PlotSpec)
if(nargin < 2)
  PlotSpec = 'cell';
end

XPad = 0.075;
YPad = 0.075;

if(StringCheck(PlotSpec, 'cell'))
  SpecType = 0;
elseif(StringCheck(PlotSpec, 'model'))
  SpecType = 1;
else
  error(['Can''t plot spectrum for: ', PlotSpec])
end

%First get some info about the map
MinFreq = Inf;
MaxFreq = 15;
for m = 1:length(Map.g_syn)
  if(SpecType == 0)
    Spec = Map.Results(m).CellReal.SlowWave.Spectrum;
  else
    Spec = Map.Results(m).ModelSlow.Spectrum;
  end
  TempMinFreq = min(Spec.Freq);
    
  if(TempMinFreq < MinFreq)
    MinFreq = TempMinFreq;
  end
end

MinPower = Inf;
MaxPower = -Inf;
for m = 1:length(Map.g_syn)
  if(SpecType == 0)
    Spec = Map.Results(m).CellReal.SlowWave.Spectrum;
  else
    Spec = Map.Results(m).ModelSlow.Spectrum;
  end
  [PFact, PInd] = GetPFact(Spec, MinFreq, MaxFreq);
  TempMinPower = min(Spec.Power(PInd)) / PFact;
  TempMaxPower = max(Spec.Power(PInd)) / PFact;
    
  if(TempMinPower < MinPower)
    MinPower = TempMinPower;
  end
  if(TempMaxPower > MaxPower)
    MaxPower = TempMaxPower;
  end
end


g_synList = unique(Map.g_syn);
g_hList = unique(Map.g_h);
NumX = length(g_synList);
NumY = length(g_hList);
StartX = g_synList(1);
DeltaX = g_synList(2) - StartX;
MaxX = DeltaX * NumX;
StartY = g_hList(1);
DeltaY = g_hList(2) - StartY;
MaxY = DeltaY * NumY;

Width = (1.0 - XPad) / NumX;
Height = (1.0 - YPad) / NumY;

yTop = sqrt(MaxPower) * 1.05;
Title = sprintf('Spectrum map for %s %s', Map.Type, PlotSpec);
h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');
clf;
hold on
h_axes = axes('Position', [0 0 1 1], 'Visible', 'off');
for m = 1:length(Map.g_syn)
  x = XPad + (1 - XPad) * (Map.g_syn(m) - StartX) / MaxX;
  y = YPad + (1 - YPad) * (Map.g_h(m) - StartY) / MaxY;
  switch(Map.Cats(m))
   case 0, C_m = 'k';
   case 1, C_m = 'g';
   case 2, C_m = 'b';
   case 3, C_m = 'r';
   otherwise, error('Invalid Category');
  end
  
  if(SpecType == 0)
    Cell = Map.Results(m).CellReal;
    Spec = Cell.SlowWave.Spectrum;
    if(Cell.Burst.Freq > 0)
      SlowFreq = Cell.Burst.Freq;
      SpikeFreq = Cell.Burst.SpikeFrequencies.Mean;
    else
      SlowFreq = Cell.SlowWave.Freq;
      SpikeFreq = Cell.Spike.Freq;
    end
  else
    ModelSlow = Map.Results(m).ModelSlow;
    Spec = ModelSlow.Spectrum;
    SlowFreq = ModelSlow.Freq;
    SpikeFreq = 0;
  end

  axes('Position', [x, y, Width, Height])
  
  [PFact, PInd] = GetPFact(Spec, MinFreq, MaxFreq);
  if(SlowFreq > 0)
    plot([SlowFreq, SlowFreq], [0, sqrt(MaxPower)], 'k--')
    hold on
  end
  if(SpikeFreq > 0)
    plot([SpikeFreq, SpikeFreq], [0, sqrt(MaxPower)], 'c:', ...
	 'LineWidth', 3)
    hold on
  end
  plot(Spec.Freq(PInd), sqrt(Spec.Power(PInd) / PFact), [C_m, '-'])
  hold off
  xlim([0, MaxFreq])
  ylim([0, yTop])
  set(gca, 'xTickLabel', {}, 'yTickLabel', {})
end
%keyboard
set(h, 'CurrentAxes', h_axes)
text(.5, 0.01,'g_s_y_n', 'FontSize', 18) 
text(0.5 * XPad, .5, 'g_h', 'Rotation', 90, 'FontSize', 18)
%xlabel('g_s_y_n')
ylabel('g_h')
hold off
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PFact, PInd] = GetPFact(Spec, MinFreq, MaxFreq)
PInd = find(Spec.Freq <= MaxFreq);
SPower = Spec.Power(PInd);
SFreq = Spec.Freq(PInd);
PFact = trapz(SFreq, SPower);
return
