%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script visualizes the energy loss cross-section as a function of
%%% indicent energy and energy loss
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

%% 1.0 Extracting the data
optData     =   scattdata.optical.inel_dcsdata;

figure(3002);
imagesc(optData.Elossmat)

eMat4EnergyDuffXSection =   optData.E'*ones([1 size(optData.Elossmat,2)]);
figure(3003);
imagesc(eMat4EnergyDuffXSection)

figure(3004);
s_eLossXSec = surf(eMat4EnergyDuffXSection,optData.Elossmat,optData.dsigdE);
xlabel('Incident Energy (eV)')
ylabel('Energy Loss (eV)')
s_eLossXSec.EdgeColor = 'none';

%% Normalize to the max of each incident energy
dsigdE_N = optData.dsigdE./...
    (max(optData.dsigdE,[],2)*ones([1 size(optData.dsigdE,2)]));
figure(3005)
s_eLossXSec_N = surf(eMat4EnergyDuffXSection,optData.Elossmat,dsigdE_N);
xlabel('Incident Energy (eV)')
ylabel('Energy Loss (eV)')
s_eLossXSec_N.EdgeColor = 'none';
colormap(jet(512))
