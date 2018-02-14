function DrawCorrelateProps(Organized, Props)

Conditions = Organized.Conditions;
ExpNums = unique(Organized.ExpNum);
PropsExpNums = cat(1, Props.ExpNum);
for n = 1:length(Conditions)
  Condition = Conditions{n};
  TypeName = sprintf('Type_%s', Condition);
  h = NamedFigure(sprintf('Correlate %s', Condition));
  set(h, 'WindowStyle', 'docked');
  x_ptx = [];
  x_mod = [];
  y = [];
  Vals = Organized.(TypeName);
  for m = 1:length(ExpNums)
    Ind = find(PropsExpNums == ExpNums(m));
    if(length(Ind) == 0)
        continue
    end
    P1 = Props(Ind(1));
    P2 = Props(Ind(2));

    if(strcmp(Condition, 'ptx') | strcmp(Condition, P1.Modulator))
        xVal_ptx = .5 * ((P1.VSpike_ptx - P1.VRest_ptx) / P1.R_ptx + ...
        	   (P2.VSpike_ptx - P2.VRest_ptx) / P2.R_ptx);
        xVal_mod = .5 * ((P1.VSpike_mod - P1.VRest_mod) / P1.R_mod + ...
        	   (P2.VSpike_mod - P2.VRest_mod) / P2.R_mod);
    else
      continue;
    end
    x_ptx = [x_ptx, xVal_ptx];
    x_mod = [x_mod, xVal_mod];
    
    Ind = find(Organized.ExpNum == ExpNums(m));
    yVal = sum(Vals(Ind) > 0);
    y = [y, yVal];
    
    %plot(xVal, yVal, 'b.');
  end
  if(strcmp(Condition, 'ptx'))
      plot(x_ptx, y, 'b.');
  else
      plot(x_mod, y, 'b.');
      hold on
      plot(x_ptx, y, 'r.');
      hold off
      legend('ptx', Condition)
  end
  xlabel('Current Needed to Excite (nA)')
  ylabel('Number of Bursters')
end

return