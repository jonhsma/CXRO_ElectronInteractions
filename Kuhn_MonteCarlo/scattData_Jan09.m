% This function initializes the scattering data in the configuration of
% Jan09 configuration of the simulator.

function scattdata = scattData_Jan09()
    % Global functions
    addpath('..\GlobalFunctions');

    % Components
    %addpath('..\CommonComponents');

    cross_sect_data_path='..\Traj_MonteCarlo\CrossSect_Data\';

    scattdata.vibr=load([cross_sect_data_path 'Vibr_Khakoo_2013_x50.mat']);
    %%% Determining the vibrational chanel of scattering
    scattdata.vibr.datasrc='Khakoo';
    %%% The Frolich IMFP function
    scattdata.vibr.imfp_func=@ephscatt;
    %%% Optical data path
    optdata_path='..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';

    scattdata.optical=load([optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData.mat']);
    pathname='..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';
    filename='DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,200].mat';
    scattdata.optical.inel_dcsdata  = load([pathname filename]);

    scattdata.E_inel_thr            = min(scattdata.optical.E);

    %% Stone wall
    % STONEWALL OPTIONS
    %%%% The energy below which the electron would activate an acid and die
    scattdata.stoneWall.CUTOFF      =   5;
    %%%% The imfp of the stone wall. Should be small if active
    scattdata.stoneWall.IMFP        =   0.001;  
    %%%% The reaction radius of the stone wall
    scattdata.stoneWall.ACID_REACTION_RADIUS   =   3;  

end