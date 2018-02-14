function varargout = FitChiSquared(fHandle, StartParams, ParamRanges, ...
				   xData, yData, ...
				   Tol, fTol, DerivTol, ...
				   Verbose, varargin)
% [Params, ChiSquared, NumEval] = FitChiSquared(fHandle, StartParams, ...
%                                               ParamRanges, ...
%				                xData, yData, ...
%				                Tol, fTol, DerivTol, ...
%				                Verbose, varargin)
% Performs a chi-squared fit of the funtion specified by fHandle,
% to xData and yData.
%   INPUT:
%     -fHandle:      Fitting function.  If the function returns
%        derivatives they should be the 2nd returned argument (and
%        they should only be returned if nargout == 2)
%     -StartParams:  Starting parameters
%     -ParamRanges:  NumParams x 2 matrix of ranges (min, max) for
%        parameters
%     -xData:        Points to evaluate fHandle
%     -yData:        Data to fit
%     -Tol:          Tolerance on parameter values
%     -fTol:         Tolerance of ChiSquared
%     -DerivTol:     Step size for finite-difference derivatives
%     -Verbose:      Set to 1 for verbose minimization, 0 otherwise
%     -varargin:     Extra parameters will be passed to fHandle
%   OUTPUT:
%     -Params:       Values of fit parameters
%     -ChiSquared:   Chi^2 per degree of freedom
%     -NumEval:      Number of function evaluations

if(nargout > 3)
  error('Error, too many input arguments.');
end

if(size(xData, 1) > 1)
  xData = xData';
end
if(size(yData, 1) > 1)
  yData = yData';
end

ChiSqrd = @ErrorFunc;
[Params, ChiSquared, NumEval] = DownHill(ChiSqrd, StartParams, ParamRanges, ...
					 Tol, fTol, DerivTol, ...
					 Verbose, xData, yData, fHandle, ...
					 varargin{:});

switch(nargout)
 case 1
   varargout = {Params};
 case 2
   varargout = {Params, ChiSquared};
 case 3
   varargout = {Params, ChiSquared, NumEval};
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Err, ErrGrad] = ErrorFunc(Params, xData, yData, fHandle, varargin)
DegOfFreedom = length(xData) - length(Params);
if(nargout == 1)
  Mod = fHandle(Params, xData, varargin{:});
  if(size(Mod, 1) > 1)
    Mod = Mod';
  end
  if(size(yData, 1) > 1)
    yData = yData';
  end
  Mod = Mod - yData;
  Err = Mod * Mod' / DegOfFreedom;
else
  [Mod, fGrad] = fHandle(Params, xData, varargin{:});
  if(size(Mod, 1) > 1)
    Mod = Mod';
  end
  if(size(yData, 1) > 1)
    yData = yData';
  end
  if(size(fGrad, 1) ~= size(Mod, 2))
    fGrad = fGrad';
  end
  
  
  Mod = Mod - yData;
  Err = Mod * Mod' / DegOfFreedom;
  ErrGrad = 2 * (Mod * fGrad) / DegOfFreedom;
end
return