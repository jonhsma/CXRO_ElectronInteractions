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

scattdata.optical.inel_dcsdata  = load([pathname filename]);

%% The reference scattering cross-section
addiParam = [];
addiParam.onlyimfp =0;

object = genrandEloss_OptData(scattdata.optical,80,scattdata.optical.inel_dcsdata,addiParam);

%% The cross-section data
energy = 80;

dcs_data =  scattdata.optical.inel_dcsdata;
idx      =  find(dcs_data.E(1,:)==energy);

if ~isempty(idx)
    thetavec=dcs_data.thetamat(idx,:); % The theta axis
    dsigdOmega=dcs_data.dsigdOmega(idx,:);
else
    [~,idx]      =  min(abs(dcs_data.E(1,:)-energy));
    thetavec=dcs_data.thetamat(idx,:); % The theta axis
    dsigdOmega=dcs_data.dsigdOmega(idx,:);
end
figure(1000)
hold off
plot(thetavec,dsigdOmega/max(dsigdOmega));
hold on
plot(thetavec,sin(thetavec).*dsigdOmega/max(sin(thetavec).*dsigdOmega),'.');
title('Cross-section per unit solid angle as a function of \theta');
legend('d\sigma/d\Omega','Corresponding sine-weighted P(\theta)');
xlabel('\theta angle')   

%% The Monte Carlo--Toss the Dice
nTrials =   1000;

%%% Initialize everything to increase speed
theta_S   =   zeros([1 nTrials]); 
phi_S     =   zeros([1 nTrials]);
resVec_S  =   zeros([3 nTrials]);

theta_W   =   zeros([1 nTrials]); 
phi_W     =   zeros([1 nTrials]);
resVec_W  =   zeros([3 nTrials]);


genLoopTimer = tic;


calibrationTest = 0;

fprintf('The dice throw begins\n');
% Parallelization slows it down
for ii = 1:nTrials
    if calibrationTest
        theta_S = acos(2*rand-1);
        phi_S   = 2*rand*2*pi; 
    end
    object  = genrandEloss_OptData(scattdata.optical,energy,scattdata.optical.inel_dcsdata,addiParam);
    theta_S   =   object.theta;
    phi_S     =   object.phi;
    
    resVec_S(:,ii) =...
        [sin(theta_S)*cos(phi_S);...
        sin(theta_S)*sin(phi_S);...
        cos(theta_S)];
end
%%% The second run with the sine modified PDF
for ii = 1:nTrials
    if calibrationTest
        theta_W = acos(2*rand-1);
        phi_W   = 2*rand*2*pi; 
    end
    object  = genrandEloss_OptData_JHM(scattdata.optical,energy,scattdata.optical.inel_dcsdata,addiParam);
    theta_W   =   object.theta;
    phi_W     =   object.phi;
    
    resVec_W(:,ii) =...
        [sin(theta_W)*cos(phi_W);...
        sin(theta_W)*sin(phi_W);...
        cos(theta_W)];
end
toc(genLoopTimer)
%% Sampling
phiRes = pi/32;
phiQuery = 0:phiRes:2*pi;
thetaRes = pi/32;
thetaQuery = 0:thetaRes:pi;

dotProdLowLimit = 0.999;

thetaResults_S = zeros([1 length(thetaQuery)]);
thetaResults_W = zeros([1 length(thetaQuery)]);
finalResults_S = zeros([length(phiQuery) length(thetaQuery)]);
finalResults_W = zeros([length(phiQuery) length(thetaQuery)]);
dotProductArray = zeros([1 nTrials]);


counter = 1;
fprintf('Now the sampling bit.\n')
samplingTimer = tic;
% Parallelization slows it down
for phi = phiQuery
    for tt = 1:length(thetaQuery) %%Preperation for parallelization
        refVec =...
            [sin(thetaQuery(tt))*cos(phi);...
            sin(thetaQuery(tt))*sin(phi);...
            cos(thetaQuery(tt))];
        dotProductArray = refVec'*resVec_S;
        thetaResults_S(tt) = length(find(dotProductArray>dotProdLowLimit));
    end
    finalResults_S(counter,:) = thetaResults_S;
    counter = counter+1;
end
counter = 1;
fprintf('Now the sampling bit.\n')
% Parallelization slows it down
for phi = phiQuery
    for tt = 1:length(thetaQuery) %%Preperation for parallelization
        refVec =...
            [sin(thetaQuery(tt))*cos(phi);...
            sin(thetaQuery(tt))*sin(phi);...
            cos(thetaQuery(tt))];
        dotProductArray = refVec'*resVec_W;
        thetaResults_W(tt) = length(find(dotProductArray>dotProdLowLimit));
    end
    finalResults_W(counter,:) = thetaResults_W;
    counter = counter+1;
end
toc(samplingTimer)
%% Display the result
finalResults_S = finalResults_S/sum(sum(finalResults_S));
finalResults_W = finalResults_W/sum(sum(finalResults_W));
figure(1001);
hold off;
surf(thetaQuery,phiQuery,finalResults_S);
axis([0,pi/2,0,2*pi,0,max(max(finalResults_S))]);
caxis([0 max(max(finalResults_S))]);


figure(1002);
hold off
semilogy(thetavec,dsigdOmega/max(dsigdOmega));
hold on;
semilogy(thetaQuery,mean(finalResults_S,1)/max(mean(finalResults_S,1)));
semilogy(thetaQuery,mean(finalResults_W,1)/max(mean(finalResults_W,1)));
xlabel('\theta angle');
ylabel('Relative scattering cross-section');
legend('Data','MC-Orginal PD','MC-Sine-Weighted PD')
        
        
%% This is a not-so-automated script to demonstrate the effect of MC sampling
figure(7004);
hold on
semilogy(thetaQuery,mean(finalResults_S,1)/max(mean(finalResults_S,1)))
legend('Data','Large Sampling Window','Small Sampling Window')



