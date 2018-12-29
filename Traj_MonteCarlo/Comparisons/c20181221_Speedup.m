%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script compares the output before and after the speed up (which
%%% end up around 1.5x, instead of the 3x suggested by the serial timed
%%% run) took place Dec20, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureWindowStyle','Docked')
%% Loading all energies of the pre-speed up data
loadSlow.runsToCompare =...
    {'20181217_SEM_LEonly';'20181217_SEM_LEonly';...
    '20181217_SEM_LEonly';'20181217_SEM_LEonly';...
    '20181217_SEM_LEonly';'20181217_SEM_LEonly'};
loadSlow.incEnergy =...
    {'30','40','50','65','80','85'};
loadSlow.doseStr  = {'0.00';'0.00';'0.00';'0.00';'0.00';'0.00'};
loadSlow.nTrials  = {88,88,88,88,88,88};
loadSlow.res      =   0.5;
loadSlow.cuOption =   'cdf';
loadSlow.distOption =   'pdf';
loadSlow.rangeLimit =   8;
%% Actually loading the slow data
slowResult = Comparison_ScatteringSim_Acid(loadSlow);
%% Redraw
plotAcidDistCurves(slowResult,loadSlow)

%% Loading all energies of the post-speed up data
loadFast.runsToCompare =...
    {'20181220_SEM_LE_PostSpeedUp';'20181220_SEM_LE_PostSpeedUp';...
    '20181220_SEM_LE_PostSpeedUp';'20181220_SEM_LE_PostSpeedUp';...
    '20181220_SEM_LE_PostSpeedUp';'20181220_SEM_LE_PostSpeedUp'};
loadFast.incEnergy =...
    {'30','40','50','65','80','85'};
loadFast.doseStr  = {'0.00';'0.00';'0.00';'0.00';'0.00';'0.00'};
loadFast.nTrials  = {88,88,88,88,88,88};
loadFast.res      =   0.5;
loadFast.cuOption =   'cdf';
loadFast.distOption =   'pdf';
loadFast.rangeLimit =   8;
%% Actually loading the plasmon data
fastResult = Comparison_ScatteringSim_Acid(loadFast);
%% Redraw
plotAcidDistCurves(slowResult,loadFast)

%% 85 eV comparison
plotAcidDistCurves_Z({slowResult{6},fastResult{6}})

%% 65 eV comparison
plotAcidDistCurves_Z({slowResult{4},fastResult{4}})

%% Does the difference matter?
comConfig.res = 1;
comConfig.incEnergy = {'85','85','65','65','40','40','30','30'};
comConfig.nTrials = {88,88,88,88,88,88,88,88};
comConfig.runsToCompare = {'85eV-slow','85eV-fast',...
    '65eV-slow','65eV-fast'...
    '40eV-slow','40eV-fast'...
    '30eV-slow','30eV-fast'};
plotAcidDistCurves_Z({slowResult{6},fastResult{6},...
    slowResult{4},fastResult{4},...
    slowResult{2},fastResult{2},...
    slowResult{1},fastResult{1},},comConfig)