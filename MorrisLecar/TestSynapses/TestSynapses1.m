function TestSynapses(FileName, g_syn, g_h)
if(nargin < 3)
  FileName = '/mnt/dwidget/758_91/758_91_LP_0056.abf';
  g_syn = 70;
  g_h = 25;
end
[t, VReal, VModel, I_Injected] = GetData(FileName);

SimData = SimDynamicClamp(t, VReal, VModel, g_syn, g_h, 0);

VModelCalc = round(2 * SimData.VModelCalc) / 2;
I_InjectedCalc = round(16 * SimData.I_InjectedCalc) / 16;

t = 0.001 * t;
NamedFigure('ML Voltage')
hold off
plot(t, VModel, 'b.')
hold on
plot(t, VModelCalc, 'r.')
hold off

NamedFigure('Injected Current')
hold off
plot(t, I_Injected, 'b.')
hold on
plot(t, I_InjectedCalc, 'r.')
hold off

Err = (VModelCalc - VModel);% ./ (VModel + 1e-3);
InjErr = (I_InjectedCalc - I_Injected);
NamedFigure('ML Voltage/Current Error')
hold off
plot(t, Err, 'b.')
hold on
plot(t, InjErr, 'g.') 
hold off

CumulativeErrVals = sort(abs(Err));
NamedFigure('ML Error Cumulative Distribution')
plot(CumulativeErrVals, linspace(1, 0, length(Err)))
fprintf('Min error magnitude: %g\n', CumulativeErrVals(1))
fprintf('Median error magnitude: %g\n', median(abs(Err)))
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