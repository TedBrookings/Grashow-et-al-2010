function CompareXDistYDist(x, y, label)
numCells = size(x, 1);
numCompare = numCells * (numCells - 1) / 2;
dx = zeros(numCompare, 1);
dy = zeros(numCompare, 1);

m = 1;
for n1 = 1:numCells
  for n2 = (n1+1):numCells
    dx(m) = dist(x(n1,:), x(n2,:));
    dy(m) = dist(y(n1,:), y(n2,:));
    m = m + 1;
  end
end
[coef, pVal] = Pearson(dx, dy);
fprintf('%s: (PearsonCorrelation, pVal) = (%g, %g)\n', label, coef, pVal)

range = [min(min(dx), min(dy)), max(max(dx), max(dy))];

titleStr = sprintf('IP vs Map Distances: %s', label);
h = NamedFigure(titleStr);
set(h, 'WindowStyle', 'docked')
hold off
plot(range, range, 'r-')
hold on
plot(dx, dy, 'b.')
hold off

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = dist(a, b)
delta = a - b;
d = sqrt(delta * (delta'));
return