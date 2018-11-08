%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file keeps track of what comparisons I make in the process of fine
%%% tuning the low energy scatteirng

%% Nov 8, 2018 80 eV. Before/After introducing low energy events
%%% This is already after I turned off the parts of the code where an acid is
%%% generated wherever a scattering event occurs.
cnfg_Nov8_80eV.runsToCompare =...
    {'LowEnergyEnabled_20181107';...
    'Ref_20181103'};
cnfg_Nov8_80eV.incEnergy =...
    {'80.00';'80.00'};
cnfg_Nov8_80eV.doseStr  = {'0.00';'0.00'};
cnfg_Nov8_80eV.nTrials  = {500,1000};
cnfg_Nov8_80eV.res      =   0.5;
cnfg_Nov8_80eV.cuOption =   'cumcount';
cnfg_Nov8_80eV.distOption =   'count';
cnfg_Nov8_80eV.rangeLimit =   20;

Comparison_ScatteringSim_Acid(cnfg_Nov8_80eV);

%% Nov 8, 2018 30 eV. Before/After introducing low energy events
%%% This is already after I turned off the parts of the code where an acid is
%%% generated wherever a scattering event occurs.
cnfg_Nov8_30eV.runsToCompare =...
    {'LowEnergyEnabled_20181107';...
    'Ref_20181103'};
cnfg_Nov8_30eV.incEnergy =...
    {'30.00';'30.00'};
cnfg_Nov8_30eV.doseStr  = {'0.00';'0.00'};
cnfg_Nov8_30eV.nTrials  = {500,500};
cnfg_Nov8_30eV.res      =   0.5;
cnfg_Nov8_30eV.cuOption =   'cumcount';
cnfg_Nov8_30eV.distOption =   'count';
cnfg_Nov8_30eV.rangeLimit =   20;

Comparison_ScatteringSim_Acid(cnfg_Nov8_30eV);