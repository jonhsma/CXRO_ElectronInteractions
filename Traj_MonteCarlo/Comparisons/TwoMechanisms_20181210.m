%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script compares the two runs with two activation mechanisms. The
% first one has only plasmon scattering and the second one has only low
% energy mechanisms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureWindowStyle','Docked')
%% Loading all energies of the plasmon data 
loadPlasmon.runsToCompare =...
    {'20181207_LowEnergyOff';'20181207_LowEnergyOff';...
    '20181207_LowEnergyOff';'20181207_LowEnergyOff';...
    '20181207_LowEnergyOff';'20181207_LowEnergyOff'};
loadPlasmon.incEnergy =...
    {'20','30','40','50','65','80'};
loadPlasmon.doseStr  = {'0.00';'0.00';'0.00';'0.00';'0.00';'0.00'};
loadPlasmon.nTrials  = {1100,1100,1100,1100,1100,1100};
loadPlasmon.res      =   0.5;
loadPlasmon.cuOption =   'cdf';
loadPlasmon.distOption =   'pdf';
loadPlasmon.rangeLimit =   5;
%% Actually loading the plasmon data
plasmonResult = Comparison_ScatteringSim_Acid(loadPlasmon);
%% Redraw
plotAcidDistCurves(plasmonResult,loadPlasmon)
%% Some additional display
addpath('components\')
mOption.runsToCompare =...
    {'20181207_LowEnergyOff';'20181207_LowEnergyOff';...
    '20181207_LowEnergyOff'};
mOption.incEnergy =...
    {'20','50','80'};
mOption.doseStr  = {'0.00';'0.00';'0.00'};
mOption.nTrials  = {1100,1100,1100};
mOption.res      =   0.5;
mOption.cuOption =   'cdf';
mOption.distOption =   'pdf';
mOption.rangeLimit =   8;
plotAcidDistributions({plasmonResult{1},...
    plasmonResult{4},plasmonResult{6}},mOption)
%%  Efficiency and spread
plotEfficiency(plasmonResult,loadPlasmon)
plotSpread(plasmonResult,loadPlasmon)

%% Loading all energies of the low energy data
loadStoneWall.runsToCompare =...
    {'20181207_LowEnergyOnly';'20181207_LowEnergyOnly';...
    '20181207_LowEnergyOnly';'20181207_LowEnergyOnly';...
    '20181207_LowEnergyOnly';'20181207_LowEnergyOnly'};
loadStoneWall.incEnergy =...
    {'20','30','40','50','65','80'};
loadStoneWall.doseStr  = {'0.00';'0.00';'0.00';'0.00';'0.00';'0.00'};
loadStoneWall.nTrials  = {1100,1100,1100,1100,1100,1100};
loadStoneWall.res      =   0.5;
loadStoneWall.cuOption =   'cdf';
loadStoneWall.distOption =   'pdf';
loadStoneWall.rangeLimit =   8;
%% Actually loading the low energy data
stoneWallResult = Comparison_ScatteringSim_Acid(loadStoneWall);
%% Redraw
plotAcidDistCurves(stoneWallResult,loadStoneWall)
%% Some additional display
addpath('components\')
mOption.runsToCompare =...
    {'20181207_LowEnergyOnly';'20181207_LowEnergyOnly'};%...
    %'20181207_LowEnergyOnly';'20181207_LowEnergyOnly'};
mOption.incEnergy =...
    {'30','80'};
mOption.doseStr  = {'0.00';'0.00';'0.00';'0.00'};
mOption.nTrials  = {1100,1100,1100,1100};
mOption.res      =   0.5;
mOption.cuOption =   'cdf';
mOption.distOption =   'pdf';
mOption.rangeLimit =   10;
plotAcidDistCurves({stoneWallResult{2},...
    stoneWallResult{6}},mOption)
%%  Efficiency and spread
plotEfficiency(stoneWallResult,loadStoneWall)
plotSpread(stoneWallResult,loadStoneWall)
%% Model
plot(xx,min(floor(xx/20),(xx+11)/31)./xx*2)
plot(xx,min(ceil(xx/20),(xx+11)/31)./xx*2)
legend('Plasmon','LE Attachment','Model')


%% 20181212 For Mentor Presentation
%% 30 eV comparison
addpath('components\')
mOption.runsToCompare =...
    {'20181207_LowEnergyOff';'20181207_LowEnergyOnly';...
    '20181207_LowEnergyOff'};
mOption.incEnergy =...
    {'30','30'};
mOption.doseStr  = {'0.00';'0.00'};
mOption.nTrials  = {1100,1100};
mOption.res      =   0.5;
mOption.cuOption =   'cdf';
mOption.distOption =   'pdf';
mOption.rangeLimit =   8;
mOption.legendArray =   {'PlasmonLoss';'LowEnergyAttachment'};
plotAcidDistributions({plasmonResult{2},stoneWallResult{2}},mOption)
%% 80 eV comparison
addpath('components\')
mOption.runsToCompare =...
    {'20181207_LowEnergyOff';'20181207_LowEnergyOnly';...
    '20181207_LowEnergyOff'};
mOption.incEnergy =...
    {'80','80'};
mOption.doseStr  = {'0.00';'0.00'};
mOption.nTrials  = {1100,1100};
mOption.res      =   0.5;
mOption.cuOption =   'cdf';
mOption.distOption =   'pdf';
mOption.rangeLimit =   8;
mOption.legendArray =   {'PlasmonLoss';'LowEnergyAttachment'};
plotAcidDistributions({plasmonResult{6},stoneWallResult{6}},mOption)