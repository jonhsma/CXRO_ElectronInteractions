% Graphics for advanced litho paper
%% load data
load('eScan_20190201')
corr_MFP = exp(smooth(eScan_20190201.Lcoef(1,:),30));
%% Run KhunLength.mlx to prepare the workspace
%% Extract correlation adjusted mean free path
binEdges = 4:4:92;
binCenter = (binEdges(1:end-1)+binEdges(2:end))/2;
fontSize = 16;

figure(4001)
hold off
hin = histogram([data_09.f_09_up_5.Ein],binEdges);
hold on;
hpe = histogram([data_09.f_09_escape.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts');
title('Photoelectron Energy Spectra')
legend('Internal','Emission');

figure(7200)
subplot(1,2,2)
hold off
plot(binCenter,hin.Values./hpe.Values,'b-','LineWidth',2)
hold on
plot(binCenter,1.5*hin.Values./hpe.Values.*interp1(energyScale,imfp,binCenter),...
    'r-','LineWidth',2)
plot(binCenter,0.5*hin.Values./hpe.Values.*interp1(eScale,corr_MFP,binCenter),...
    'm-','LineWidth',2)
xlabel('KE(eV)','FontSize',fontSize);
ylabel('Ratio','FontSize',fontSize);
title({'Proportionality'},'FontSize',fontSize)
set(gca,'FontSize', fontSize)
legend({'Internal/Emission',...
    '1.5\times Internal/(Emission/MFP)',...
    '0.5\times Internal/(Emission/Corr-Adj-MFP)'},'FontSize',fontSize*0.75);
axis([-inf inf 0 inf])
%% Displaying the correlation adjusted mean-free path
%% Plot the mean free path
figure(7200)
subplot(1,2,1)
hold off
fontSize = 16;
plot(eScale,corr_MFP,'y-','LineWidth',2);
hold on;
plot(energyScale, imfp,'g-','LineWidth',2)
plot(energyScale, imfp_opt_a,'-.')
plot(energyScale, imfp_vibr_a,'-.')
plot(energyScale, imfp_stnw_a,'r-.')
legend('Correlation-Adjusted MFP','Total MFP','Plasmon','Vibrational','Low energy absorption');
axis([-inf 90 0 3.5])
set(gca,'FontSize',fontSize)
xlabel('KE (eV)')
ylabel('Mean free path (nm)');
