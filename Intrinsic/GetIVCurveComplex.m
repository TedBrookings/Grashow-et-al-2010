function IV = GetIVCurve(t, v, I)
%First find when the current turns on and off
ABS_I = abs(I);
[Max1, Ind1] = max(ABS_I);
[MaxI, Ind2] = max(Max1);

Ind = find(ABS_I(:,Ind2) > .1 * MaxI);

On = Ind(1) + 1;
Off = Ind(end) - 1;

%Get the resting membrane potential
V0 = mean(median(v(2:(On-2),:)));

%Loop through the traces
NumTraces = size(v, 2);
ShowPlot = true;
NoShape = true;
for n = 1:NumTraces
  %Check to see if it's spiking
  Analyze = AnalyzeWaveform2(t(On:Off), v(On:Off,n), ShowPlot, NoShape);
  
  if(length(Analyze.Spike.Times) > 0)
    %If it's spiking, we're done with the list
    break
  end
  %Otherwise, get the current, voltage, and time constants
  [IV_I(n), IV_V(n), DeltaV1(n), Tau1(n), DeltaV2(n), Tau2(n)] ...
      = GetIVAndTau(t, I, v, V0, On, Off, n);
end

figure
plot(t, v)
keyboard

LineVals = polyfit(IV_I, IV_V, 1);

IV.I = IV_I;
IV.V = IV_V;
IV.DeltaV1 = DeltaV1;
IV.Tau1 = Tau1;
IV.DeltaV2 = DeltaV2;
IV.Tau2 = Tau2;
IV.R = LineVals(1);
IV.VRest = V0;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = GetIVAndTau(t, I, v, V0, On, Off, Trace)
IV_I = mean(I(On:Off,Trace));

xData = t(On:Off);
yData = v(On:Off, Trace);

DeltaV = yData(end) - V0;
if(DeltaV < 0)
  VRange1 = [1.5 * DeltaV, .7 * DeltaV];
  VRange2 = [.9 * DeltaV, -.5 * DeltaV];
else
  VRange1 = [.7 * DeltaV, 1.5 * DeltaV];
  VRange2 = [-.5 * DeltaV, .9 * DeltaV];
end

fHandle = @VCurve;
StartParams = [V0, .8 * DeltaV, .005, .2 * DeltaV, .1];
ParamRanges = [[V0, V0]; VRange1; [.001, .3]; VRange2; [.02, 2]];

DeltaT = t(2) - t(1);
Tau1 = DeltaT / log(abs((yData(3) - yData(1))/(yData(2) - yData(1))));
DeltaV1 = yData(2) - yData(1);
DeltaV2 = DeltaV - DeltaV1;
Tau2 = 10 * Tau1;
%StartParams = [V0, DeltaV1, Tau1, DeltaV2, Tau2];

VRange1 = sort([-DeltaV * .5, DeltaV * 1.5]);
VRange2 = sort([-DeltaV * .5, DeltaV * 1.5]);
TauRange1 = [Tau1 / 10, Tau1 * 10];
TauRange2 = [Tau2 / 10, Tau2 * 10];
%ParamRanges = [[V0, V0]; VRange1; TauRange1; VRange2; TauRange2];


Tol = 1e-4;
fTol = 1e-5;
DerivTol = 1e-7;

Verbose = false;

try
  Params = FitChiSquared(fHandle, StartParams, ParamRanges, ...
			 xData, yData, ...
			 Tol, fTol, DerivTol, ...
			 Verbose);
catch
  
  disp('Error fitting in GetIVCurve.m, GetIVAndTau()')
  keyboard
end
DebugFit = true;
if(DebugFit)
  Params
  hold off
  plot(xData, yData)
  hold on
  yFit = VCurve(Params, xData);
  plot(xData, yFit, 'r-');
  
  keyboard
end

  
DeltaV1 = Params(2);
Tau1 = Params(3);
DeltaV2 = Params(4);
Tau2 = Params(5);

IV_V = V0 + DeltaV1 + DeltaV2;

varargout = {IV_I, IV_V, DeltaV1, Tau1, DeltaV2, Tau2};
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V, varargout] = VCurve(Parameters, t)
%Parameters(1) = V0  (fixed by data)
%Parameters(2) = DeltaV1
%Parameters(3) = Tau1
%Parameters(4) = DeltaV2
%Parameters(5) = Tau2

Temp1 = 1 - exp(-(t - t(1)) / Parameters(3));
Temp2 = 1 - exp(-(t - t(1)) / Parameters(5));
V = Parameters(1) + Parameters(2) * Temp1 ...
    + Parameters(4) * Temp2;
return
if(nargout > 1)
  Deriv = [ones(size(t)), Temp1, ...
	   Parameters(2) / Parameters(3)^2 * (Temp1 .* (t - t(1))), ...
	   Temp2, ...
	   Parameters(4) / Parameters(5)^2 * (Temp2 .* (t - t(1)))];
  varargout = {Deriv};
end
return