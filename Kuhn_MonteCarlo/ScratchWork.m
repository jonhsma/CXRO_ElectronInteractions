% This is a script for the Monte carlo that's built to find the Kuhn IMFP
%% This is a scratch file so DON'T BATCH RUN ANYTHING FROM IT

%% Paths
addpath('../Traj_MonteCarlo/Discrete_Energy_Losses_Approach_2')
addpath('../Traj_MonteCarlo/Discrete_Energy_Losses_Approach_2/StoneWall')
addpath('../Traj_MonteCarlo/Discrete_Energy_Losses_Approach_2/OptDataScatt')
addpath('../Traj_MonteCarlo/Discrete_Energy_Losses_Approach_2/VibrationScatt')

%% Scattering data
scattdata = scattData_Jan09;

%% Try to run one trajectory
initConfig.energy = 92;
initConfig.xyz = [0 0 0];
initConfig.theta_in = 0; % initial condition for theta in z
initConfig.phi_in = 0; % initial condition for the polar angle
initConfig.nSteps = 5000;

nSteps          =   initConfig.nSteps;
results_singlePass = TrajectoryFollowerConstEnergy(initConfig,...
    scattdata,...
    12/120*6.02*1e23,...
    @arrayAllocator);
%% Test the correlation anaysis function
tic
[testSummary,testM] = trajCorr(results_singlePass);
toc
% it took 0.1 second so parallelism is not needed.
%% How did that go
figure(7300)
hold off
semilogx(testSummary.mean)
errorbar(testSummary.mean,(testSummary.std)./sqrt(testSummary.N'+1));
set(gca,'XScale','Log')
xlabel('Separation (n,events)')
ylabel('Correlation');
title('Correlation <cos\theta_{i,i+n}> as a function of separation')
grid on
figure(7301)
hold off
semilogx(testSummary.mean)
errorbar(testSummary.mean,(testSummary.std)./sqrt(testSummary.N'+1)*sqrt(20));
set(gca,'XScale','Log')
xlabel('Separation (n,events)')
ylabel('Correlation');
title({'Correlation <cos\theta_{i,i+n}>';' Error enlarged by 20^{1/2}'})
grid on
%% Test the length analysis function
[testSummary_L,testM_L] = trajLengthAna(results_singlePass);
%% The results
error = (testSummary_L.ini.std)./sqrt(testSummary_L.ini.N+1);
range = 1:5000;
figure(7201)
hold off
loglog(range,testSummary_L.ini.mean(range))
errorbar(range,testSummary_L.ini.mean(range),...
    error(range));
hold on
loglog(range,sqrt(testSummary_L.ini.ms(range)))
set(gca,'YScale','Log')
set(gca,'XScale','Log')
title('<R>,<R^2>^{1/2} vs N')
xlabel('Number of steps');
ylabel('nm')
legend('<R>','<R^2>^{1/2}')
grid on

figure(7203)
hold off
range = 30:400;
plot(range,testSummary_L.ini.mean(range))
hold on
lm = fitlm(log(range),log(testSummary_L.ini.mean(range)),...
    'Weights',testSummary_L.ini.mean(range)./error(range));% Error has to be transformed as y is logged
plot(range,exp(lm.predict))
axis([0 max(range) -inf inf])
title('<R> vs N and fitting')
xlabel('Number of steps');
ylabel('<R>')
legend('Sim',strcat(num2str(exp(lm.Coefficients.Estimate(1))),'N^{',...
    num2str(lm.Coefficients.Estimate(2)),'}'));

%% Test the energy scan thingy
esTest = energyScan([80 50 30],5000,scattdata);

%% Full run on Feb 01 12:50 am
eScale = 5:2:90;
eScan_20190131 = energyScan(eScale,5000,scattdata);

%% Analyze that set of data
[X,Y] = meshgrid(5:2:90,1:5000);
figure(7210)
surf(X,Y,eScan_20190131.Lmean)
shading interp
xlabel('KE(eV)');ylabel('N');zlabel('<R> (nm)');

% A zoom in
limit = 300;
figure(7211)
[X,Y] = meshgrid(5:2:90,1:limit);
surf(X,Y,eScan_20190131.Lmean(1:limit,:))
shading interp
xlabel('KE(eV)');ylabel('N');zlabel('<R> (nm)');

% The Kuhn mean free path
figure(7212)
hold off
errorbar(eScale,exp(eScan_20190131.Lcoef(1,:)),...
    exp(eScan_20190131.Lcoef(1,:)).*eScan_20190131.Lcoef_se(2,:))
xlabel('KE(eV)')
ylabel('\lambda_{kuhn}(nm)')
title('Kuhn Inelastic Mean Free Path')
hold on
errorbar(eScale,eScan_20190131.Lcoef(2,:),eScan_20190131.Lcoef_se(2,:))
% Smooth the K-IMFP and append it
kimfp_s = smooth(exp(eScan_20190131.Lcoef(1,:)),7);
plot(eScale,kimfp_s)
legend('MonteCarlo Results','\gamma in <R> = \lambda N^{\gamma}','Smoothed \lambda(E)')

%% This part requires running KuhnLeng.mlx
figure(7213)
hold off;
hin = histogram([data_09.f_09_up_5.Ein],binEdges);
hold on;
hpe = histogram([data_09.f_09_escape.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts');
title('Photoelectron Energy Spectra')
legend('Internal','Emission');

figure(7214)
hold off
plot(binCenter,hin.Values./hpe.Values)
hold on
% The huristic version
plot(binCenter,0.5*hin.Values./hpe.Values.*interp1(energyScale,kuhn_IMFP,binCenter))
% The montecarlo version
plot(binCenter,0.5*hin.Values./hpe.Values.*interp1(eScale,kimfp_s,binCenter))
xlabel('KE(eV)');
ylabel('Counts');
title({'Proportionality'})
legend('Internal/Emission','Interanl/Emission(Huristic-IMFP corrected) * 1/2',...
    'Interanl/Emission(MC Kuhn-IMFP corrected)');
axis([-inf inf 0 inf])

%% Summary
% The Kuhn mean free path
figure(7212)
hold off
errorbar(eScale,exp(eScan_20190131.Lcoef(1,:)),...
    exp(eScan_20190131.Lcoef(1,:)).*eScan_20190131.Lcoef_se(2,:))
xlabel('KE(eV)')
ylabel('\lambda_{kuhn}(nm)')
title('Kuhn Inelastic Mean Free Path')
hold on
errorbar(eScale,eScan_20190131.Lcoef(2,:),eScan_20190131.Lcoef_se(2,:))
% Smooth the K-IMFP and append it
kimfp_s = smooth(exp(eScan_20190131.Lcoef(1,:)),7);
plot(eScale,kimfp_s)
range = 3:90;
plot(energyScale(range),kuhn_IMFP(range))
legend('MonteCarlo Results','\gamma in <R> = \lambda N^{\gamma}',...
    'Smoothed \lambda(E)','Huristic \lambda_{Kuhn}')

%% Investigate the 'toll' of takeing mean instead of taking RMS
% Redo the summary with more data
% But test it first. Test the energy scan thingy
esTest = energyScan([92 50 30],5000,scattdata);

%% Full run on Feb 01 6:03 pm
eScale = 5:0.5:92; %176 energies = 8 iterations
eScan_20190201 = energyScan(eScale,5000,scattdata);
%% Visualization
range = 1:300;
[X,Y] = meshgrid(eScale,range);
figure(7210)
surf(X,Y,sqrt(eScan_20190201.Lms(range,:)));
shading interp
set(gca,'YScale','Log')
set(gca,'ZScale','Log')
%% INFP comparison
figure(7213)
hold off
plot(eScale,smooth(eScan_20190201.Lcoef_ms(2,:),10).^0.5)
hold on
plot(eScale,exp(smooth(eScan_20190201.Lcoef(1,:),5)))
ratio = smooth(eScan_20190201.Lcoef_ms(2,:),10).^0.5./...
    exp(smooth(eScan_20190201.Lcoef(1,:),5));
plot(eScale,ratio)
% The mean square fitting always give larger IMFP
mean(eScan_20190201.Lcoef_ms(1,:))
% The intercept is negative
%% Try forcing the intercept zero
range = 30:400;
for ei = 1:length(eScale)
        error = eScan_20190201.Lsdse(:,ei)./sqrt(length(eScan_20190201.Lmean)+1);
        lm = fitlm(log(range),log(eScan_20190201.Lms(range,ei)),...
            'Weights',eScan_20190201.Lms(range,ei)./error(range));%,...
            %'Intercept',false);
    eScan_20190201_aux.coef(:,ei) = lm.Coefficients.Estimate;
    eScan_20190201_aux.coef_se(:,ei) = lm.Coefficients.SE;
end
%% INFP comparison
boxcar = 20;
figure(7213)
hold off
plot(eScale,smooth(eScan_20190201.Lcoef_ms(2,:),boxcar).^0.5)
hold on
plot(eScale,exp(smooth(eScan_20190201.Lcoef(1,:),boxcar)))
plot(eScale,exp(smooth(eScan_20190201_aux.coef(1,:),boxcar)/2))
ratio1 = smooth(eScan_20190201.Lcoef_ms(2,:),boxcar).^0.5./...
    exp(smooth(eScan_20190201.Lcoef(1,:),boxcar));
plot(eScale,ratio1)
legend('Linear fit <R^2>','Log fit <R>','Log fit <R^2>','\lambda_{lin}/\lambda_{log}')
ylabel('nm');
xlabel('KE(eV)');
title({'Comparison on fitting methods';...
    strcat('Curves smoothed by a ',num2str(boxcar),' boxcar')});



