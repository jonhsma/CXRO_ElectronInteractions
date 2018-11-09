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

%% Nov 8, 2018 80 eV. Before/After fixing low energy
%%% This is already after I increased the vibrational cross-section and
%%% turned off "free energy" acid generation events
cnfg_Nov8_80eV2.runsToCompare =...
    {'LowEnergyEnabled_20181107';...
    'LowEnergyFixed_20181108'};
cnfg_Nov8_80eV2.incEnergy =...
    {'80.00';'80.00'};
cnfg_Nov8_80eV2.doseStr  = {'0.00';'0.00'};
cnfg_Nov8_80eV2.nTrials  = {500,500};
cnfg_Nov8_80eV2.res      =   0.5;
cnfg_Nov8_80eV2.cuOption =   'cumcount';
cnfg_Nov8_80eV2.distOption =   'count';
cnfg_Nov8_80eV2.rangeLimit =   20;

Comparison_ScatteringSim_Acid(cnfg_Nov8_80eV2);

%% Nov 8, 2018 80 eV. Before enabling LE/After fixing low energy
%%% This is already after I increased the vibrational cross-section and
%%% turned off "free energy" acid generation events
cnfg_Nov8_80eV3.runsToCompare =...
    {'Ref_20181103';...
    'LowEnergyFixed_20181108'};
cnfg_Nov8_80eV3.incEnergy =...
    {'80.00';'80.00'};
cnfg_Nov8_80eV3.doseStr  = {'0.00';'0.00'};
cnfg_Nov8_80eV3.nTrials  = {500,500};
cnfg_Nov8_80eV3.res      =   0.5;
cnfg_Nov8_80eV3.cuOption =   'cumcount';
cnfg_Nov8_80eV3.distOption =   'pdf';
cnfg_Nov8_80eV3.rangeLimit =   10;

Comparison_ScatteringSim_Acid(cnfg_Nov8_80eV3);
%% Nov 8, 2018 30 eV. Before enabling LE/After fixing low energy
%%% This is already after I increased the vibrational cross-section and
%%% turned off "free energy" acid generation events
cnfg_Nov8_30eV3.runsToCompare =...
    {'Ref_20181103';...
    'LowEnergyFixed_20181108'};
cnfg_Nov8_30eV3.incEnergy =...
    {'30.00';'30.00'};
cnfg_Nov8_30eV3.doseStr  = {'0.00';'0.00'};
cnfg_Nov8_30eV3.nTrials  = {500,500};
cnfg_Nov8_30eV3.res      =   0.5;
cnfg_Nov8_30eV3.cuOption =   'cumcount';
cnfg_Nov8_30eV3.distOption =   'pdf';
cnfg_Nov8_30eV3.rangeLimit =   10;

Comparison_ScatteringSim_Acid(cnfg_Nov8_30eV3);
%% Nov 8, 2018 30 eV. Before enabling LE/After fixing low energy
%%% This is already after I increased the vibrational cross-section and
%%% turned off "free energy" acid generation events
cnfg_Nov8.runsToCompare =...
    {'LowEnergyFixed_20181108';...
    'LowEnergyFixed_20181108'};
cnfg_Nov8.incEnergy =...
    {'80.00';'30.00'};
cnfg_Nov8.doseStr  = {'0.00';'0.00'};
cnfg_Nov8.nTrials  = {500,500};
cnfg_Nov8.res      =   0.5;
cnfg_Nov8.cuOption =   'cumcount';
cnfg_Nov8.distOption =   'pdf';
cnfg_Nov8.rangeLimit =   10;

Comparison_ScatteringSim_Acid(cnfg_Nov8);