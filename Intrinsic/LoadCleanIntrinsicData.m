function [t, v, I, Okay] = LoadCleanIntrinsicData(FileName, CellName)
%[t, v, I, Okay] = LoadCleanIntrinsicData(FileName, CellName)
%  INPUT PARAMETERS:
%   -FileName:  name with full path to .abf file
%  OUTPUT PARAMETERS:
%   -t is time trace in ms
%   -v is voltage tracein mV
%   -I is current trace in nA
%   -Okay is a boolean that is true if there are no errors, false otherwise

if(nargin ~= 2)
  error(['Incorrect number of input arguments.  ', ...
	 'Run "help LoadCleanIntrinsicData"'])
end

AbfS = LoadAbf(FileName);
t = AbfS.Time * 1000;  %convert to ms

FNames = fieldnames(AbfS.Units);
Current = [];
Voltage = [];
for n = 1:length(FNames)
  Unit = AbfS.Units.(FNames{n});
  if(strcmp(Unit, 'mV'));
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
  if(~StringCheck(FNames{Voltage}, CellName))
    Okay = false;
    return
  elseif(StringCheck(FNames{Current}, CellName))
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
  if(StringCheck(FNames{Voltage(1)}, CellName))
    v_n = 1;
    VOther_n = 2;
  else
    v_n = 2;
    VOther_n = 1;
  end
  if(StringCheck(FNames{Current(1)}, CellName))
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
    VOther = AbfS.Data.(FNames{Voltage(VOther_n)});
    if(std(VOther) > 2)
      %Bridge mode, we requested the "bad" cell
      Okay = false;
      return
    else
      %DCC
      Okay = true;
      DCC = true;
    end
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

if(DCC)  %interpolate to DCC sampling frequency
  [t, v, I, DCC_Info] = CleanDCC(t, v, I);
  if(~isfinite(DCC_Info.DCC_Freq))
    disp('Warning:  weird DCC signal!')
    DCC_Info
  else
    disp(sprintf('DCC Freq = %g kHz', DCC_Info.DCC_Freq))
  end
end

%Low-pass filter with a characteristic time ~0.5 ms
%SmoothTime = 0.5;
SmoothTime = 1.0;
if(t(2) - t(1) < SmoothTime)
  n = 2;
  while(t(n) - t(1) < SmoothTime)
    n = n + 2;
  end
  [B,A] = butter(2, 2 / n,'low');
  v = filtfilt(B, A, v);
end
return