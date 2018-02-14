function PlotNarrow()
%Loading this isn't necessary for now:
%{
ExperimentsExists = evalin('base', 'exist(''ML_Experiments'')');
if(ExperimentsExists)
  Experiments = evalin('base', 'ML_Experiments');
else
  Experiments = LoadAllExperiments('MorrisLecarFolders.txt');
  assignin('base', 'ML_Experiments', Experiments);
end
%}
close all

%Draw burst fraction maps:  Change to plot narrow fraction instead...
%DrawBurstFraction(Experiments, 'ML');

%Load intrinsic properties
Intrinsics = LoadIPxls;

[narrowProps, narrowLabels] = GetNarrowNetworks;

%Select which kinds of plots to make:
PlotIndividual = false;
PlotImportance = false;

%Choose which ScoreHandle to use:
ScoreHandle = @ZScore;
%ScoreHandle = @RankScore;

%Choose the number of shuffled trials (should be large or zero):
NumShuffledTrials = 000;

fitOptions.FitPooled = false;
fitOptions.FitSubTypes = true;
fitOptions.IPvsIP = false;
fitOptions.CompareDist = true;
fitOptions.FitNetProps = false;

%Do the analysis with the selected map properties
ResultStruct = CorrelateIntrinsics(Intrinsics, ...
				   narrowProps, narrowLabels, ...
				   ScoreHandle, fitOptions);
if(NumShuffledTrials > 0)
  ResultStruct.ShufflePVals = zeros(ResultStruct.NumFits, ...
				    NumShuffledTrials);
  ResultStruct.ShuffleRSquared = zeros(ResultStruct.NumFits, ...
				       NumShuffledTrials);  
  fprintf('\nRunning shuffled trials: ')
  for n = 1:NumShuffledTrials
    ResultStruct.ShuffleNum = n;
    ResultStruct = CorrelateIntrinsics(Intrinsics, ...
				       narrowProps, narrowLabels, ...
				       ScoreHandle, fitOptions, ...
				       ResultStruct);

    if(n == 21 || (n < 100 && mod(n - 20, 30) == 1) || ...
       (n > 100 && mod(n, 20) == 1))
      fprintf('\n %g', n)
    else
      fprintf(' %g', n)
    end
  end
  fprintf('\n')
else
  ResultStruct.ShufflePVals = [];
  ResultStruct.ShuffleRSquared = [];
end

CorrectAndDisplay(ResultStruct, PlotIndividual, PlotImportance);
return
