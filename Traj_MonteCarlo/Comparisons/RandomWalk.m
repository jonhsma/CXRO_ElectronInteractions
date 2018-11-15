%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Charaterizing the random walk runs to make sure that the gears are
%%% working fine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5 and 10 eV comparison, With low energy events
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