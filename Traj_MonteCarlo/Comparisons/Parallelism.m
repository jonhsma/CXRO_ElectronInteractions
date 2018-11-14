%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parallelism is implemented around Nov 10 2018.
%%% The purpose of the script is to investigate whether the parallilization
%%% is causing unwanted changes or abnormalities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 80 eV comparison, No Low energy events
cnfg_Nov13.runsToCompare =...
    {'Ref_20181103';...
    'Parallel_22W_Ref_20181109'};
cnfg_Nov13.incEnergy =...
    {'80.00';'80.00'};
cnfg_Nov13.doseStr  = {'0.00';'0.00'};
cnfg_Nov13.nTrials  = {1000,1000};
cnfg_Nov13.res      =   0.5;
cnfg_Nov13.cuOption =   'cumcount';
cnfg_Nov13.distOption =   'pdf';
cnfg_Nov13.rangeLimit =   10;

Comparison_ScatteringSim_Acid(cnfg_Nov13);

%% 30 eV comparison, No Low energy events
cnfg_Nov13.runsToCompare =...
    {'Ref_20181103';...
    'Parallel_22W_Ref_20181109'};
cnfg_Nov13.incEnergy =...
    {'30.00';'30.00'};
cnfg_Nov13.doseStr  = {'0.00';'0.00'};
cnfg_Nov13.nTrials  = {1000,1000};
cnfg_Nov13.res      =   0.5;
cnfg_Nov13.cuOption =   'cumcount';
cnfg_Nov13.distOption =   'pdf';
cnfg_Nov13.rangeLimit =   10;

Comparison_ScatteringSim_Acid(cnfg_Nov13);

%% 80 eV comparison, With low energy events
cnfg_Nov13.runsToCompare =...
    {'LowEnergyFixed_20181108';...
    'Parallel_22W_LowE_80_Ref_20181113'};
cnfg_Nov13.incEnergy =...
    {'80.00';'80.00'};
cnfg_Nov13.doseStr  = {'0.00';'0.00'};
cnfg_Nov13.nTrials  = {500,500};
cnfg_Nov13.res      =   0.5;
cnfg_Nov13.cuOption =   'cumcount';
cnfg_Nov13.distOption =   'pdf';
cnfg_Nov13.rangeLimit =   10;

Comparison_ScatteringSim_Acid(cnfg_Nov13);

%% 30 eV comparison, With low energy events
cnfg_Nov13.runsToCompare =...
    {'LowEnergyFixed_20181108';...
    'Parallel_22W_LowE_30_Ref_20181113'};
cnfg_Nov13.incEnergy =...
    {'30.00';'30.00'};
cnfg_Nov13.doseStr  = {'0.00';'0.00'};
cnfg_Nov13.nTrials  = {500,500};
cnfg_Nov13.res      =   0.5;
cnfg_Nov13.cuOption =   'cumcount';
cnfg_Nov13.distOption =   'pdf';
cnfg_Nov13.rangeLimit =   10;

Comparison_ScatteringSim_Acid(cnfg_Nov13);
