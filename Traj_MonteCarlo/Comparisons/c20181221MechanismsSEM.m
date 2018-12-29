%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script compares the two activation mechanisms in a vertical
% incidence context. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureWindowStyle','Docked')
%% Loading all energies of the Low energy  data
loadLE.runsToCompare =...
    {'20181217_SEM_LEonly';'20181217_SEM_LEonly';...
    '20181217_SEM_LEonly';'20181217_SEM_LEonly';...
    '20181217_SEM_LEonly';'20181217_SEM_LEonly'};
loadLE.incEnergy =...
    {'30','40','50','65','80','85'};
loadLE.doseStr  = {'0.00';'0.00';'0.00';'0.00';'0.00';'0.00'};
loadLE.nTrials  = {88,88,88,88,88,88};
loadLE.res      =   0.5;
loadLE.cuOption =   'cdf';
loadLE.distOption =   'count';
loadLE.rangeLimit =   5;
%% Actually loading the LE data
leResult = Comparison_ScatteringSim_Acid(loadLE);
%% Redraw
plotAcidDistCurves_Z(leResult,loadLE)

%% Loading all energies of the Low energy  data
loadHE.runsToCompare =...
    {'20181217_SEM_LEoff';'20181217_SEM_LEoff';...
    '20181217_SEM_LEoff';'20181217_SEM_LEoff';...
    '20181217_SEM_LEoff';'20181217_SEM_LEoff'};
loadHE.incEnergy =...
    {'30';'40';'50';'65';'80';'85'};
loadHE.doseStr  = {'0.00';'0.00';'0.00';'0.00';'0.00';'0.00'};
loadHE.nTrials  = {88,88,88,88,88,88};
loadHE.res      =   0.5;
loadHE.cuOption =   'cdf';
loadHE.distOption =   'count';
loadHE.rangeLimit =   5;
%% Actually loading the LE data
heResult = Comparison_ScatteringSim_Acid(loadHE);
%% Redraw
plotAcidDistCurves_Z(heResult,loadHE)

%% Comparison
comConfig.runsToCompare =...
    {'LE-80';'HE-80';'LE-40';'HE-40'};
comConfig.incEnergy =...
    {'80';'80';'40';'40'};
comConfig.doseStr  = {'0.00';'0.00';'0.00';'0.00'};
comConfig.nTrials  = {88,88,88,88};
comConfig.res      =   0.5;
comConfig.cuOption =   'cdf';
comConfig.distOption =   'count';
comConfig.rangeLimit =   5;
plotAcidDistCurves_Z({leResult{6},heResult{6},...
    leResult{2},heResult{2}},comConfig);

%% Efficiency Comparison
plotEfficiency(heResult,loadHE,25)
plotEfficiency(leResult,loadLE,25)


