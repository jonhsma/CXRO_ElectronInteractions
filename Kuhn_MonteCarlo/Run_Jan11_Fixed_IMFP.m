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
eScan_fixedIMFP = energyScan(eScale,5000,scattdata,@TrajectoryFollowerConstEnergyConstIMFP);
%% Result analysis
boxcar = 20;
% The mean displacement fit. The prefactor is log(lambda)
lambda_mean_r_log = exp(eScan_fixedIMFP.Lcoef(1,:));
lambda_mean_r_log_s = smooth(lambda_mean_r_log,boxcar);
% The mean displacement^2 fit. The coefficient for x is lambda^2
lambda_mean_rsq_linear = sqrt(eScan_fixedIMFP.Lcoef_ms(2,:));
lambda_mean_rsq_linear_s = smooth(lambda_mean_rsq_linear,boxcar);
figure(3000)
hold off
plot(eScale,lambda_mean_r_log_s);
hold on
plot(eScale,lambda_mean_rsq_linear_s);
xlabel('KE (eV)')
ylabel('Correlation adjusted mean free path')
legend({'Log fit','Linear fit'},'Location','southeast')
%% Does it work
% Load Jan 11 data manually
binEdges    = 4:4:90;
fontSize    =   16;
binCenters  = (binEdges(1:end-1)+binEdges(2:end))/2;
figure(3001)
hold off
hpe = histogram([escapeEvents_11.Ein],binEdges);
hold on
hin = histogram([f_11_up_5.Ein],binEdges);

rawRatio = hin.Values./hpe.Values;
correctedRatio = hin.Values./hpe.Values.*interp1(eScale,lambda_mean_r_log_s,binCenters);

figure(14400)
hold off
plot(binCenters,rawRatio,'LineWidth',2);
hold on
plot(binCenters,0.6*correctedRatio,'LineWidth',2);
xlabel('KE (eV)','FontSize',fontSize)
ylabel('Ratio','FontSize',fontSize);
title({'Results for Correlation adjusted mean free path (\lambda_{cor})method',...
    'For Jan 11 data where a fixed IMFP is used',...
    '\lambda_{cor} is obtained by MC'},...
    'FontSize',fontSize);
axis([5 85 0 3.5])
legend({'Internal/Emission','0.6\times Internal/(Emission/\lambda_{cor})'},...
    'Location','northeast',...
    'FontSize',fontSize)
set(gca,'FontSize',fontSize)