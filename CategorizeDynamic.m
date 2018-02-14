function [Cat, CatString] = CategorizeDynamic(Analysis)
%[Cat, CatString] = CategorizeDynamic(Analysis)
% categorizes results from AnalyzeDynamic and puts into one of four
% mutually-exclusive categories:
%  0 - Silent
%  1 - Assymmetrical
%  2 - Spiking
%  3 - Half-Center

PowCutoff = 2;
if(Analysis.Cell0.Spike.Freq <= 0 & Analysis.Cell1.Spike.Freq <= 0)
  Cat = 0;
  CatString = 'Silent';
elseif(Analysis.Cell0.Spike.Freq <= 0 | Analysis.Cell1.Spike.Freq <= 0)
  Cat = 1;
  CatString = 'Assymmetrical';
elseif(Analysis.HalfCenter.Freq == 0 | Analysis.HalfCenter.ExclusionFact < .1)
  Cat = 2;
  CatString = 'Spiking';
else
  Cat = 3;
  CatString = 'Half-Center';
end
return