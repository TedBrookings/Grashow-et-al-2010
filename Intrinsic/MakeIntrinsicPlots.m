function MakeIntrinsicPlots(IVList, RampList, FIList)

MakeComparePlots(IVList, 'IV', {'R'}, 'Input Impedence');
MakeComparePlots(IVList, 'IV', {'VRest'}, 'VRest');
MakeComparePlots(IVList, 'IV', {'VIntercept'}, 'Projected VRest');
MakeComparePlots(RampList, 'Ramp', {'PreMaxCurve', 'V'}, 'Spike Threshold');
%MakeComparePlots(FIList, 'FI', {'SpikeHeight', 'H'}, 'Spike Height');
%MakeComparePlots(FIList, 'FI', {'SpikeWidth', 'W'}, 'Spike Width');

MakeXYPlots(RampList, 'Ramp', {'PreMaxCurve', 'I'}, 'ISpike', ...
	    IVList, 'IV', {'R'}, 'R', 'ptx');
MakeXYPlots(RampList, 'Ramp', {'SpikeHeight'}, 'Spike Height', ...
	    IVList, 'IV', {'R'}, 'R', 'ptx');
MakeXYPlots(RampList, 'Ramp', {'SpikeHeight'}, 'Spike Height', ...
	    RampList, 'Ramp', {'PreMaxCurve', 'I'}, 'ISpike', 'ptx');


MakeXYPlots(RampList, 'Ramp', {'PreMaxCurve', 'I'}, 'ISpike', ...
	    IVList, 'IV', {'R'}, 'R', 'ptx');


MakeXYPlots(FIList, 'FI', {'PreMaxCurve', 'V'}, 'Voltage', ...
	    FIList, 'FI', {'PreMaxCurve', 'K'}, 'Curvature');
MakeXYPlots(FIList, 'FI', {'Deriv', 'V'}, 'Voltage', ...
	    FIList, 'FI', {'Deriv', 'D'}, 'Max. Derivative');
MakeXYPlots(FIList, 'FI', {'SpikeHeight', 'V'}, 'Voltage', ...
	    FIList, 'FI', {'SpikeHeight', 'H'}, 'Spike Height');
MakeXYPlots(FIList, 'FI', {'SpikeWidth', 'V'}, 'Voltage', ...
	    FIList, 'FI', {'SpikeWidth', 'W'}, 'Spike Width');

MakeXYPlots(FIList, 'FI', {'I'}, 'Current', ...
	    FIList, 'FI', {'F'}, 'Freq');
MakeComparePlots(FIList, 'FI', {'Slope'}, 'FI Slope');
MakeComparePlots(FIList, 'FI', {'Offset'}, 'FI Offset');
MakeLineComparePlots(FIList)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeComparePlots(ExpList, LookIn, Subfields, TitleStr)

TypeList = fieldnames(ExpList(1).Files);
%First make the pooled histogram plots
for Type = 1:length(TypeList)
  Vals = [];
  TypeStr = TypeList{Type};
  for Num = 1:length(ExpList)
    Vals = [Vals, GetVal(ExpList(Num), LookIn, TypeStr, Subfields)];
  end
  Ind = strfind(TypeStr, '_') + 1;
  FullTitle = sprintf('(%s) %s', TypeStr(Ind:end), TitleStr);
  h = NamedFigure(FullTitle);
  set(h, 'WindowStyle', 'docked');
  hist(Vals, 20);
end

%Now make the comparison plots
for Type = 1:length(TypeList)
  TypeStr = TypeList{Type};
  if(strcmp(TypeStr, 'file_ptx'))
    continue
  end
  Val_ptx = [];
  Val_mod = [];
  for Num = 1:length(ExpList)
    Temp_ptx = GetVal(ExpList(Num), LookIn, 'file_ptx', Subfields);
    if(length(Temp_ptx) == 0)
      continue;
    end
    Temp_mod = GetVal(ExpList(Num), LookIn, TypeStr, Subfields);
    if(length(Temp_mod) == 0)
      continue;
    end
    Val_ptx = [Val_ptx, Temp_ptx];
    Val_mod = [Val_mod, Temp_mod];
  end
  Ind = strfind(TypeStr, '_') + 1;
  FullTitle = sprintf('(ptx->%s) %s', TypeStr(Ind:end), TitleStr);
  h = NamedFigure(FullTitle);
  set(h, 'WindowStyle', 'docked');
  ValRange = [min([min(Val_ptx), min(Val_mod)]), ...
	      max([max(Val_ptx), max(Val_mod)])];
  plot(ValRange, ValRange, 'r-');
  hold on
  plot(Val_ptx, Val_mod, 'bo');
  hold off
  xlim(ValRange);
  ylim(ValRange);
  axis square
  xlabel('ptx', 'FontSize', 16)
  ylabel(TypeStr(Ind:end), 'FontSize', 16)
  title(TitleStr, 'FontSize', 16)
  
  [Sig, pVal] = ttest(Val_mod - Val_ptx);
  disp(sprintf('T-test for %s:  pVal = %g', FullTitle, pVal))
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Val = GetVal(ExpStruct, LookIn, TypeStr, Subfields)
if(length(ExpStruct.Files.(TypeStr)) == 0)
  Val = [];
  return
end
Val = ExpStruct.(LookIn).(TypeStr).(Subfields{1});
for n = 2:length(Subfields)
  Val = Val.(Subfields{n});
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeXYPlots(XList, XLookIn, XSubfields, XTitle, ...
		     YList, YLookIn, YSubfields, YTitle, ValidType)

TitleStr = sprintf('%s vs. %s', XTitle, YTitle);
TypeList = fieldnames(XList(1).Files);

if(nargin == 9)
  switch(ValidType)
   case 'ptx', MarkerStyle = 'bo';
   case '5ht', MarkerStyle = 'rx';
   case 'oxo', MarkerStyle = 'g+';
   otherwise,
    error(sprintf('Invalid type:  %s', ValidType))
  end
else
  MakeXYPlots(XList, XLookIn, XSubfields, XTitle, ...
	      YList, YLookIn, YSubfields, YTitle, ...
	      'ptx');
  hold on
  MakeXYPlots(XList, XLookIn, XSubfields, XTitle, ...
	      YList, YLookIn, YSubfields, YTitle, ...
	      '5ht');
  MakeXYPlots(XList, XLookIn, XSubfields, XTitle, ...
	      YList, YLookIn, YSubfields, YTitle, ...
	      'oxo');
  legend('ptx', '5ht', 'oxo');
  hold off
  return
end

for Type = 1:length(TypeList)
  TypeStr = TypeList{Type};
  Ind = strfind(TypeStr, '_') + 1;
  if(nargin == 9)
    if(~strcmp(TypeStr(Ind:end), ValidType))
      continue;
    end
  end
  XVals = [];
  YVals = [];
  for XNum = 1:length(XList)
    X = GetVal(XList(XNum), XLookIn, TypeStr, XSubfields);
    if(length(X) == 0)
      continue;
    end
    YNum = find(strcmp({YList.ID}, XList(XNum).ID));
    if(length(YNum) == 0)
      continue;
    end
    Y = GetVal(YList(YNum), YLookIn, TypeStr, YSubfields);
    if(length(Y) == 0)
      continue;
    end
    XVals = [XVals, X];
    YVals = [YVals, Y];
  end
  %  FullTitle = sprintf('(%s) %s', TypeStr(Ind:end), TitleStr);
  FullTitle = TitleStr;
  h = NamedFigure(FullTitle);
  set(h, 'WindowStyle', 'docked');
  plot(XVals, YVals, MarkerStyle);
  xlabel(XTitle);
  ylabel(YTitle);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeLineComparePlots(ExpList)
TitleStr = 'FI behavior';

TypeList = fieldnames(ExpList(1).Files);
LookIn = 'FI';
SlopeSubfield = {'Slope'};
OffsetSubfield = {'Offset'};

%Now make the comparison plots
for Type = 1:length(TypeList)
  TypeStr = TypeList{Type};
  if(strcmp(TypeStr, 'file_ptx'))
    continue
  elseif(strcmp(TypeStr, 'file_oxo'))
    TypeColor = 'g-';
  else
    TypeColor = 'r-';
  end
  
  Slope_ptx = [];
  Offset_ptx = [];
  Slope_mod = [];
  Offset_mod = [];
  for Num = 1:length(ExpList)
    TempSlope_ptx = GetVal(ExpList(Num), LookIn, 'file_ptx', SlopeSubfield);
    if(length(TempSlope_ptx) == 0)
      continue;
    end
    TempOffset_ptx = GetVal(ExpList(Num), LookIn, 'file_ptx', ...
			    OffsetSubfield);
    if(length(TempOffset_ptx) == 0)
      continue;
    end
    
    TempSlope_mod = GetVal(ExpList(Num), LookIn, TypeStr, SlopeSubfield);
    if(length(TempSlope_mod) == 0)
      continue;
    end
    TempOffset_mod = GetVal(ExpList(Num), LookIn, TypeStr, ...
			    OffsetSubfield);
    if(length(TempOffset_mod) == 0)
      continue;
    end
    
    Slope_ptx = [Slope_ptx, TempSlope_ptx];
    Offset_ptx = [Offset_ptx, TempOffset_ptx];
    Slope_mod = [Slope_mod, TempSlope_mod];
    Offset_mod = [Offset_mod, TempOffset_mod];
  end
  Ind = strfind(TypeStr, '_') + 1;
  FullTitle = sprintf('(ptx->%s) %s', TypeStr(Ind:end), TitleStr);
  h = NamedFigure(FullTitle);
  set(h, 'WindowStyle', 'docked');

  Ind = find(isfinite(Slope_ptx) & isfinite(Offset_ptx) & ...
	     isfinite(Slope_mod) & isfinite(Offset_mod));
  clf
  hold on
  I = 0:.01:3;
  for n = Ind
    V = Slope_ptx(n) * I + Offset_ptx(n);
    V(find(V < 0)) = 0;
    plot(I, V, 'b-');
    V = Slope_mod(n) * I + Offset_mod(n);
    V(find(V < 0)) = 0;
    plot(I, V, TypeColor);
  end
  hold off
end
return