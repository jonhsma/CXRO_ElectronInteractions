%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script visualizes the IMFP as a function of energy using the
%%% scattering engines used in the scattering sim.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0.0 Path and other basic initialization
clear;
clear global
clc;
close all;
set(0,'DefaultFigureWindowStyle','docked')

% Global functions
addpath('..\..\GlobalFunctions');

cross_sect_data_path='..\..\Traj_MonteCarlo\CrossSect_Data\';

scattdata.vibr=load([cross_sect_data_path 'VibrExcit_Data_Khakoo_2013.mat']);
%%% Determining the vibrational chanel of scattering
scattdata.vibr.datasrc='Khakoo';
%%% The Frolich IMFP function
scattdata.vibr.imfp_func=@ephscatt;
%%% Optical data path
optdata_path='..\..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';

scattdata.optical=load([optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData.mat']);

%TruongData=load('TruongData_PMMA_PS.mat');

pathname='..\..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';
%{
filename='DDCSdata_Ef=0p5_Elossmin=0.001eV_Erange=[5,1000].mat';

filename='DDCSdata_withICSData_Ef=10_Elossmin=0.001eV_Erange=[5,1000].mat';

filename='DDCSdata_Fuji_Ef=15.5_Elossmin=3eV_Erange=[19,200]_EQCons=Pines.mat';
%}
filename='DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,200].mat';

scattdata.optical.inel_dcsdata  = load([pathname filename]);
scattdata.E_inel_thr            = min(scattdata.optical.E);
%% 0.1 Scattering engines
addpath('..\Discrete_Energy_Losses_Approach_2\OptDataScatt\')
addpath('..\Discrete_Energy_Losses_Approach_2\VibrationScatt\')
addpath('..\Discrete_Energy_Losses_Approach_2\StoneWall\')
addpath('..\Discrete_Energy_Losses_Approach_2\')
%% 1.0 --> Model Parameters
%%% The limit where scattering ceases
SCATTERING_LOW_ENERGY_CUTOFF    =       2;
%%% The energy where the electron enters low energy regime
%%% (Where low energy interaction is turned on)
LOW_ENERGY_BEHAVIOUR_BOUNDARY   =       20; 
%%% Low energy random walk mean free path.
LOW_ENERGY_MEAN_FREE_PATH       =       3.67;
%%% The reaction radius of PAGS
ACID_REACTION_RADIUS            =       3;
%%% Molecular density for vicrational calculations
MOLECULAR_NUMBER_DENSITY        =       1.2/120*6.02*1e23*10; % molecules/cm3
%%%% The energy below which the electron would activate an acid and die
scattdata.stoneWall.CUTOFF      =   5;
%%%% The imfp of the stone wall. Should be small if active
scattdata.stoneWall.IMFP        =   0.0001;  
%%%% The reaction radius of the stone wall
scattdata.stoneWall.ACID_REACTION_RADIUS   =   3;  

%% The IMFPs
energyScale     =   2:1:100;
imfp            =   energyScale*0;
imfp_opt_a      =   energyScale*0;
imfp_vibr_a     =   energyScale*0;
imfp_stnw_a     =   energyScale*0;

counter = 0;
for ei = energyScale   
    counter     =   counter +1;
    %% Optical IMFP
    controlparm.onlyimfp=1;
    if isfield(scattdata.optical,'inel_dcsdata')
        Elossrand_opt=genrandEloss_OptData_JHM(scattdata.optical,ei,scattdata.optical.inel_dcsdata,controlparm);
    else
        Elossrand_opt=genrandEloss_OptData_JHM(scattdata.optical,ei,controlparm);
    end
    
    if ei     >   LOW_ENERGY_BEHAVIOUR_BOUNDARY
        controlparm.onlyimfp=0;
        imfp_opt=Elossrand_opt.imfp; % imfp in nm; the above lines only calculate imfp due to the controlparms.onlyimfp line abvove.
        if isnan(imfp_opt)
            imfp_opt=Inf;
        end
    else
        imfp_opt = LOW_ENERGY_MEAN_FREE_PATH;
    end
    imfp_opt_a(counter) = imfp_opt;
    %% Vibrational IMFP
    if ei < LOW_ENERGY_BEHAVIOUR_BOUNDARY
        vibrScattCmplx     =   genrandEloss_Vibr(scattdata.vibr,ei);
        invImfp_invCM      =   vibrScattCmplx.ics*MOLECULAR_NUMBER_DENSITY;
        imfp_vibr          =   1/invImfp_invCM*1e7; % imfp in nm
    else 
        imfp_vibr          =    inf;
    end
    imfp_vibr_a(counter) = imfp_vibr;
    %% StoneWall IMFP
    stoneWallResults    = scattEngineStoneWall(ei,scattdata);
    imfp_stoneWall      =   stoneWallResults.imfp;
    imfp_stnw_a(counter) = imfp_stoneWall;
    %% Total IMFP
    imfp(counter)    =   1/(1/imfp_opt+1/imfp_vibr+1/imfp_stoneWall);
end

%% Visualization
figure(2001);
hold off;
plot(energyScale,imfp)
hold on
plot(energyScale,imfp_opt_a,'*');
plot(energyScale,imfp_vibr_a,'o');
plot(energyScale,imfp_stnw_a,'X');
title('IMFP as a function of energy')
legend('Total IMFP','Optical','Vibrational','StoneWall')