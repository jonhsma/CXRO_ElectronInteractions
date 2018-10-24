%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file tests the scattering kernel by doing a Monte Carlo
%%% simulation on input parameters and sample the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading the data
%%% Paths
optdata_path    =   '..\..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';
pathname='..\..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';

%%% Specific files
filename        =   'DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,200].mat';

scattdata.optical =...
    load([optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData.mat']);

scattdata.optical.inel_dcsdata  =load([pathname filename]);

%%
addiParam = [];
addiParam.onlyimfp =0;

object = genrandEloss_OptData(scattdata.optical,80,scattdata.optical.inel_dcsdata,addiParam);

