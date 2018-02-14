function [Cat, CatString] = CategorizeML(Analysis)
%[Cat, CatString] = CategorizeML(Analysis)
% categorizes results from AnalyzeML and puts into one of four
% mutually-exclusive categories:
%  0 - Silent
%  1 - Bio spiking, Model inhibited
%  2 - Model Wins, Bio inhibited
%  3 - Half-center

SlowWaveCutoff = 2.0;  %Minimum amplitude to be "significant"

SpikeFreq = Analysis.CellReal.Spike.Freq;
NumSpikes = length(Analysis.CellReal.Spike.Times);
if(NumSpikes < 3)
  SpikeFreq = 0;
end
BurstFreq = Analysis.CellReal.Burst.Freq;
HalfCenter = Analysis.CellReal.HalfCenter;

RealSlowWaveAmp = mean(Analysis.CellReal.SlowWave.Amplitudes);
SlowWaveCutoff = 1.0 / Analysis.CellReal.SlowWave.Corr;
if(isnan(RealSlowWaveAmp) || RealSlowWaveAmp < SlowWaveCutoff)
  RealSlowWaveAmp = 0;
end
SlowWaveCutoff = 1.0 / Analysis.ModelSlow.Corr;
ModelSlowWaveAmp = mean(Analysis.ModelSlow.Amplitudes);
if(isnan(ModelSlowWaveAmp) || ModelSlowWaveAmp < SlowWaveCutoff)
  ModelSlowWaveAmp = 0;
end

if(SpikeFreq <= 0)
  if(ModelSlowWaveAmp == 0)
    Cat = 0;
    CatString = 'Silent';
  else
    Cat = 2;
    CatString = 'Model Wins, Bio inhibited';
  end
elseif(HalfCenter)
  Cat = 3;
  CatString = 'Half-center';
else
  Cat = 1;
  CatString = 'Bio spiking, Model inhibited';
end

return