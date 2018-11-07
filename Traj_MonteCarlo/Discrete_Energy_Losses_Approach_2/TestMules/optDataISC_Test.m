%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script assess if it is faster to initialize the integral
%%% cross-sections before running the simulation than calculating then on
%%% the spot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the files
addpath('..\')
addpath('..\OptDataScatt\')

optdata_path='..\DDCSData\';
scattdata.optical=load([optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData.mat']);

pathname='..\DDCSData\';
filename='DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,200].mat';
scattdata.optical.inel_dcsdata  = load([pathname filename]);

%% The initialization function
initialization  =   tic;
optdata = optDataScatt_Init(scattdata.optical);
toc(initialization)
E = optdata.inel_dcsdata.E;
thetamat = optdata.inel_dcsdata.thetamat;
figure(3000);
s = surf(E'*ones([1 size(thetamat,2)]),thetamat,optdata.inel_dcsdata.isc);
%s.EdgeColor = 'none';