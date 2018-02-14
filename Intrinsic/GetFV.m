function GetFV(FileName);
This.WhichCell = 'top';
[T, V, I, Okay] = GetCleanData(FileName, This);
if(~Okay)
  This.WhichCell = 'bot';
  [T, V, I, Okay] = GetCleanData(FileName, This);
  if(~Okay)
    error(sprintf('Could not find data in %s.', FileName));
  end
end

PlotVar = false;
NoShape = false;
f = [];
v = [];
for n = 1:size(V,2)
  Analyze = AnalyzeWaveform2(T, V(:,n), PlotVar, NoShape);
  New_f = Analyze.Spike.Frequencies.List;
  New_v = Analyze.Spike.PreMaxCurve.V.List;
  f = [f, New_f];
  v = [v, .5 * (New_v(2:end) + New_v(1:(end-1)))];
end

myfit = polyfit(v, f, 1);
f_fit = myfit(1) * v + myfit(2);

NamedFigure('f-v plot and linear fit');
plot(v, f, 'bo');
hold on
plot(v, f_fit, 'r-');
hold off;
%title('f-v plot and linear fit');
xlabel('Voltage (mV)');
ylabel('Frequency (hZ)');


NamedFigure('Voltage Trace')
hold off
plot(T, V(:,1));
hold on
for n = 2:size(V,2)
  plot(T, V(:,n))
end
hold off
%title('Voltage Trace')
xlabel('Time (ms)')
ylabel('Voltage (mV)')

VZeroFreq = -myfit(2) / myfit(1);
Slope = myfit(1);
disp(sprintf('Extrapolated voltage of zero frequency:  %g', VZeroFreq))
disp(sprintf('Slope:  %g (Hz / mV)', Slope))
disp(sprintf('Offset:  %g (Hz)', myfit(2)))
return










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, v, I, Okay] = GetCleanData(FileName, This)
%disp(sprintf('%s, %s', FileName, This.WhichCell))
AbfS = LoadAbf(FileName);
t = AbfS.Time * 1000;  %convert to ms

FNames = fieldnames(AbfS.Units);
Current = [];
Voltage = [];
for n = 1:length(FNames)
  Unit = AbfS.Units.(FNames{n});
  if(strcmp(Unit, 'mV'))
    Voltage = [Voltage, n];
  elseif(strcmp(Unit, 'nA'))
    Current = [Current, n];
  end
end
if(length(Current) == 1 & length(Voltage == 1))
  %There are two possibilities
  %1) It's DCC, and we want the voltage that corresponds to the
  %    requested cell (top/bot).
  %2) It's two-electrode bridge mode, in which case we want the
  %    voltage trace that corresponds to the requested cell, and
  %    the current that corresponds to the other cell
  v = AbfS.Data.(FNames{Voltage});
  I = AbfS.Data.(FNames{Current});
  if(length(strfind(FNames{Voltage}, This.WhichCell)) == 0)
    Okay = false;
    return
  elseif(length(strfind(FNames{Current}, This.WhichCell)) > 0)
    Okay = true;
    DCC = true;
  else
    Okay = true;
    DCC = false;
  end
elseif(length(Current) == 2 & length(Voltage == 2))
  %There are two possibilities
  %1) It's DCC, and we want the half that corresponds to the
  %    requested cell (top/bot).
  %2) It's two-electrode bridge mode, in which case one voltage
  %    trace is good, and the other current trace corresponds to
  %    the good voltage trace; and one current trace just is a
  %    noisy zero.
  if(length(strfind(FNames{Voltage(1)}, This.WhichCell)) > 0)
    v_n = 1;
  else
    v_n = 2;
  end
  if(length(strfind(FNames{Current(1)}, This.WhichCell)) > 0)
    I_n = 1;
    IOther_n = 2;
  else
    I_n = 2;
    IOther_n = 1;
  end
  
  v = AbfS.Data.(FNames{Voltage(v_n)});
  I = AbfS.Data.(FNames{Current(I_n)});
  IOther = AbfS.Data.(FNames{Current(IOther_n)});
  if(std(IOther) < .1)
    %Bridge mode, we requested the "bad" cell
    Okay = false;
    return
  elseif(std(I) < .1)
    %Bridge mode, we requested the "good" cell
    Okay = true;
    I = IOther;
    DCC = false;
  else
    %DCC
    Okay = true;
    DCC = true;
  end
else
  error(sprintf('Problem getting data from %s', FileName))
end

if(DCC)
  [t, v, I, DCC_Info] = CleanDCC(t, v, I);
  if(~isfinite(DCC_Info.DCC_Freq))
    disp('Warning:  weird DCC signal!')
    DCC_Info
  end
end

SmoothTime = .5;
if(t(2) - t(1) < SmoothTime)
  n = 2;
  while(t(n) - t(1) < SmoothTime)
    n = n + 2;
  end
  [B,A] = butter(2, 2 / n,'low');
  v = filtfilt(B, A, v);
end
return