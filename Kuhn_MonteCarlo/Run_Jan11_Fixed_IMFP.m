%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script tests the Monte-Carlo approach in the FIXED IMFP system
% The fixed-IMFP scenario is not physical and failed the huristic
% correlated IMFP test. This script will test if the Monte Carlo approach
% works better
%% Paths
addpath('../Traj_MonteCarlo/Discrete_Energy_Losses_Approach_2')
addpath('../Traj_MonteCarlo/Discrete_Energy_Losses_Approach_2/StoneWall')
addpath('../Traj_MonteCarlo/Discrete_Energy_Losses_Approach_2/OptDataScatt')
addpath('../Traj_MonteCarlo/Discrete_Energy_Losses_Approach_2/VibrationScatt')
%% Scattering data
scattdata = scattData_Jan09;
%% The scan
eScale = 5:0.5:92; %176 energies = 8 iterations
eScan_fixedIMFP = energyScan([89 90],5000,scattdata,...
    @TrajectoryFollowerConstEnergyConstIMFP);