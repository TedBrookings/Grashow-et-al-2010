function TestSynapses(FileName, g_syn, g_h)
if(nargin < 3)
  FileName = '/mnt/dwidget/758_91/758_91_LP_0056.abf';
  g_syn = 70;
  g_h = 25;
end
[t, VReal, VModel, I_Injected] = GetData(FileName);

SimData0 = SimDynamicClamp(t, VReal, VModel, g_syn, g_h, 0);
SimData1 = SimDynamicClamp(t, VReal, VModel, g_syn, g_h, 1);
SimData2 = SimDynamicClamp(t, VReal, VModel, g_syn, g_h, 2);

VModelCalc0 = round(2 * SimData0.VModelCalc) / 2;
I_InjectedCalc0 = round(16 * SimData0.I_InjectedCalc) / 16;
VModelCalc1 = round(2 * SimData1.VModelCalc) / 2;
I_InjectedCalc1 = round(16 * SimData1.I_InjectedCalc) / 16;
VModelCalc2 = round(2 * SimData2.VModelCalc) / 2;
I_InjectedCalc2 = round(16 * SimData2.I_InjectedCalc) / 16;

t = 0.001 * t;
NamedFigure('ML Voltage')
hold off
plot(t, VModelCalc0, 'b.')
hold on
plot(t, VModelCalc1, 'r.')
plot(t, VModelCalc2, 'g.')
hold off

NamedFigure('Injected Current')
hold off
plot(t, I_InjectedCalc0, 'b.')
hold on
plot(t, I_InjectedCalc1, 'r.')
plot(t, I_InjectedCalc2, 'g.')
hold off

Err1 = (VModelCalc1 - VModelCalc0);
Err2 = (VModelCalc2 - VModelCalc0);
InjErr1 = (I_InjectedCalc1 - I_InjectedCalc0);
InjErr2 = (I_InjectedCalc2 - I_InjectedCalc0);
NamedFigure('ML Voltage Error')
hold off
plot(t, Err1, 'r.')
hold on
plot(t, Err2, 'g.') 
hold off

NamedFigure('ML Current Error')
hold off
plot(t, InjErr1, 'r.')
hold on
plot(t, InjErr2, 'g.') 
hold off

%{
CumulativeErrVals = sort(abs(Err));
NamedFigure('ML Error Cumulative Distribution')
plot(CumulativeErrVals, linspace(1, 0, length(Err)))
fprintf('Min error magnitude: %g\n', CumulativeErrVals(1))
fprintf('Median error magnitude: %g\n', median(abs(Err)))
%}
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, VReal, VModel, I_Injected] = GetData(FileName)
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
if(length(Voltage) ~= 2)
  for n = 1:length(FNames)
    disp(FNames)
  end
  error(['Incorrect number of voltage traces in ', FileName])
end

if(length(Current) ~= 1)
  for n = 1:length(FNames)
    disp(FNames)
  end
  error(['Incorrect number of current traces in ', FileName])
end

I_Injected = AbfS.Data.(FNames{Current(1)});

if(StringCheck(FNames{Voltage(1)}, 'model'))
  if(StringCheck(FNames{Voltage(2)}, 'model'))
    disp(FNames)
    error(['Two model traces found in ', FileName])
  end
  VReal = AbfS.Data.(FNames{Voltage(2)});
  VModel = AbfS.Data.(FNames{Voltage(1)});
elseif(StringCheck(FNames{Voltage(2)}, 'model'))
    VReal = AbfS.Data.(FNames{Voltage(1)});
    VModel = AbfS.Data.(FNames{Voltage(2)});
else
  disp(FNames)
  error(['No model traces found in ', FileName])
end

%[t, VReal, VModel] = CleanAndSmooth(t, VReal, VModel);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, VReal, VModel] = CleanAndSmooth(t, VReal, VModel)
%First remove DCC noise by interpolating to the DCC frequency
[t, VReal, VModel, DCC_Info] = CleanDCC(t, VReal, VModel);
if(~isfinite(DCC_Info.DCC_Freq))
  disp('Warning:  weird DCC signal!')
  DCC_Info
end
disp(sprintf('DCC freq = %g kHz', DCC_Info.DCC_Freq))

%Next low-pass filter the waveform
SmoothTime = .5;  %ms
if(t(2) - t(1) < SmoothTime)
  n = 2;
  while(t(n) - t(1) < SmoothTime)
    n = n + 2;
  end
  [B,A] = butter(2, 2 / n,'low');
  VReal = filtfilt(B, A, VReal);
  VModel = filtfilt(B, A, VModel);
end
return