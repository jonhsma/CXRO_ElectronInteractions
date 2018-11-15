%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Charaterizing the random walk runs to make sure that the gears are
%%% working fine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nov 14 5 and 10 eV comparison, Stochastic pathlength, IMPF 3.67
cnfg_Nov14.runsToCompare =...
    {'RandomWalkVerfication_20181114';...
    'RandomWalkVerfication_20181114'};
cnfg_Nov14.incEnergy =...
    {'5.00';'10.00'};
cnfg_Nov14.doseStr  = {'0.00';'0.00'};
cnfg_Nov14.nTrials  = {1100,1100};
cnfg_Nov14.res      =   0.5;
cnfg_Nov14.cuOption =   'cumcount';
cnfg_Nov14.distOption =   'pdf';
cnfg_Nov14.rangeLimit =   25;

Comparison_ScatteringSim_Acid(cnfg_Nov14);

%% Nov 15 5 and 10 eV comparison, Fixed Path length 1
cnfg_Nov15.runsToCompare =...
    {'RandomWalkVerfication_fixed_20181114';...
    'RandomWalkVerfication_fixed_20181114';
    'RandomWalkVerfication_fixed_20181114';
    'RandomWalkVerfication_fixed_20181114'};
cnfg_Nov15.incEnergy =...
    {'5.00';'10.00';'30.00';'50.00'};
cnfg_Nov15.doseStr  = {'0.00';'0.00';'0.00';'0.00';};
cnfg_Nov15.nTrials  = {1100,1100,1100,1100};
cnfg_Nov15.res      =   0.25;
cnfg_Nov15.cuOption =   'cumcount';
cnfg_Nov15.distOption =   'pdf';
cnfg_Nov15.rangeLimit =   8;

Comparison_ScatteringSim_Acid(cnfg_Nov15);