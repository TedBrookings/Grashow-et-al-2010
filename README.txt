This repository is a permanant snapshot of the code used in Grashow et al 2010:
Compensation for variable intrinsic neuronal excitability by circuit-synaptic interactions.
Grashow R, Brookings T, Marder E.
J Neurosci. 2010 Jul 7;30(27):9145-56. doi: 10.1523/JNEUROSCI.0980-10.2010.

NOTES:
-This code is not maintained, and the raw data is not stored in this repository (except "IP summary for ted.txt"). The principal value is to provide ground truth for anyone seeking to understand, validate, or compare the algorithms used in the original paper.

-This directory contains a mix of files that were used and in development.  The old spike/burst analysis method is retained as AnalyzeWaveform3.m  That is the method used in Grashow et al 2010.  AnalyzeWaveform.m is a replacement method that was in development and used in subsequent work. It is not called by any methods.

To analyze one network:
  Analysis = RunAnalyze(FileName, PlotVar)      (PlotVar = 1 to plot, 0 or omitted to not plot)
To analyze an experiment:
  AnalyzeExperiment(DirName)
To analyze all experiments:
  AnalyzeAllExperiments


To make plots for one experiment:
  PlotIndividual(DirName)
To make pooled plots
  PlotPooled
