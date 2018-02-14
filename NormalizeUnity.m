function y = NormalizeUnity(x)
% y = NormalizeUnity(x)

fprintf('mean(x) = %g\n', mean(x))
fprintf('std(x)  = %g\n', std(x))
fprintf('CoV(x)  = %g\n', std(x) / abs(mean(x)))

y = (x - min(x)) / range(x);
return