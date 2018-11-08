%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The script logs all the comparsons done with the comparison script
%%% after if it being functionalized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\')

%% Nov 5, 2018 80 eV
%%%Ensuring that the refactoring is not chaning the sim results
%%% 80 eV
cnfg_Nov5_80eV.runsToCompare =...
    {'NoCoarseGrain_0_Calib';...
    'NoCoarseGrain_3_FullFunction';...
    'Ref_20181103'};
cnfg_Nov5_80eV.incEnergy =...
    {'80.00';'80.00';'80.00'};
cnfg_Nov5_80eV.doseStr = {'0.00';'0.00';'0.00'};
cnfg_Nov5_80eV.nTrials = {1000,1000,1000};
cnfg_Nov5_80eV.res      =   0.5;
cnfg_Nov5_80eV.cuOption =   'cdf';
cnfg_Nov5_80eV.distOption =   'pdf';

Comparison_ScatteringSim_Acid(cnfg_Nov5_80eV);

%% Nov 5, 2018 30 eV
%%% Ensuring that the refactoring is not chaning the sim results
%%% 30 eV
cnfg_Nov5_80eV.runsToCompare =...
    {'NoCoarseGrain_0_Calib';...
    'Ref_20181103'};
cnfg_Nov5_80eV.incEnergy =...
    {'30.00';'30.00'};
cnfg_Nov5_80eV.doseStr = {'0.00';'0.00'};
cnfg_Nov5_80eV.nTrials = {1000,1000};
cnfg_Nov5_80eV.res      =   0.5;
cnfg_Nov5_80eV.cuOption =   'cdf';
cnfg_Nov5_80eV.distOption =   'pdf';

Comparison_ScatteringSim_Acid(cnfg_Nov5_80eV);

