%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script compares the internal flux with the photoemission flux after
% fixing the meanfreepath at 0.5 nm
% RUN IT FROM DISCRETE_ENERGY_LOSSES_APPROACH_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reload the archive if I'm restarting

outputParent    =   strcat('..\\..\\..\\..\\JonathanCodeIO_CXRO\\',...
            'ElectronInteractions\\LEEMRes\\');
outputFolder    =   '20190111_PES_FixedIMFP';
outputBasePath  =   strcat(outputParent,outputFolder,'\\');
energyScanArchive = loadEnergyArchive(outputBasePath,'ScanArchive_85eV');

%% Extract the photoemission events
[escapeEvents_11,~] = extractEscapeEvents(energyScanArchive);
disp(length(escapeEvents_11));

%% Extract the internal flux at 5 m
[f_11_up_5,~] = extractZEvent(energyScanArchive,-5,1);
%

%% Plot the histograms
binEdges = 4:4:88;
figure(7200);
hold off;
histogram([escapeEvents_11.Ein],binEdges,'Normalization','pdf');
hold on;
histogram([f_11_up_5.Ein],binEdges,'Normalization','pdf');
legend('Photoemission','InternalFlux');
title('PDFs')

figure(7201);
hold off;
h0p00 = histogram([escapeEvents_11.Ein],binEdges,'Normalization','count');
hold on;
h5p00 = histogram([f_11_up_5.Ein],binEdges,'Normalization','count');
set(gca,'YScale','Log')
title('Absolute counts')
legend('Photoemission','InternalFlux');


%% Take the ratio of the two
binCenter = (binEdges(1:end-1)+binEdges(2:end))/2;
figure(7202);
plot(binCenter,h5p00.Values./h0p00.Values);

%% Energy correlate with energy loss
figure(7203);
subplot(1,2,1);
plot([f_11_up_5.Ein],[f_11_up_5.Eloss],'r.')
xlabel('KE_{ini}(eV)');
ylabel('Energy loss (eV)');
set(gca,'YScale','Linear')

subplot(1,2,2);
plot([f_11_up_5.Eout],[f_11_up_5.Ein],'r.')
xlabel('KE_f(eV)');
ylabel('KE_i(eV)');
set(gca,'YScale','Linear')

%% Same graph, segment by scattering mechanism
figure(7204);
vibr_idx = find(strcmp({f_11_up_5.scattType},'Vibrational'));
opti_idx = find(strcmp({f_11_up_5.scattType},'Optical'));
wall_idx = find(strcmp({f_11_up_5.scattType},'StoneWall'));
hold off;
plot([f_11_up_5(vibr_idx).Ein],[f_11_up_5(vibr_idx).Eloss],'r.')
hold on;
plot([f_11_up_5(opti_idx).Ein],[f_11_up_5(opti_idx).Eloss],'b.')
xlabel('KE_{ini}(eV)');
ylabel('Energy loss (eV)');
set(gca,'YScale','Log')
legend({strcat('Vibrational :', num2str(length(vibr_idx)));...
    strcat('Optical :', num2str(length(opti_idx)))});

%% Energy distribution by scattering type
binEdges = 4:4:88;
figure(7205);
hold off;
histogram([f_11_up_5(vibr_idx).Ein],binEdges,'Normalization','count');
hold on;
histogram([f_11_up_5(opti_idx).Ein],binEdges,'Normalization','count');
histogram([f_11_up_5(wall_idx).Ein],binEdges,'Normalization','count');
set(gca,'YScale','Log');
title('Scattering counts at dpeth of 5 nm for different channels')
legend('Vibrational','Optical','Absorption')

%% Energy distribution by scattering type
figure(7206);
hold off;
histogram([f_11_up_5(vibr_idx).theta],'Normalization','pdf');
hold on;
histogram([f_11_up_5(opti_idx).theta],'Normalization','pdf');
title('Scattering angle distribution at dpeth of 5 nm for different channels')
legend('Vibrational','Optical')

%% Check the last event before ejection
[~,~,preEscapeEvents_10] = extractEscapeEvents(energyScanArchive);
%% The respective (escaping) energy distribution
vibr_idx_p = find(strcmp({preEscapeEvents_10.scattType},'Vibrational'));
opti_idx_p = find(strcmp({preEscapeEvents_10.scattType},'Optical'));
prim_idx_p = find(strcmp({preEscapeEvents_10.scattType},'NewTraj'));

figure(7207);
subplot(1,2,1)
hold off;
hpv = histogram([preEscapeEvents_10(vibr_idx_p).Eout],binEdges);
hold on
hpp = histogram([preEscapeEvents_10(prim_idx_p).Eout],binEdges);
hpo = histogram([preEscapeEvents_10(opti_idx_p).Eout],binEdges);

set(gca,'YScale','Log');
title('Last scattering event of escaped electrons')
subplot(1,2,2)
hold off;
bar(binCenter',[hpv.Values',hpp.Values',hpo.Values'],'Stacked')
legend('Vibrational','Primary/Secondary','Optical');
title('Last scattering event of escaped electrons')

%% Energy spectrum of escaped electrons after removing optical events
figure(7208)
subplot(1,2,1);
hold off;
h_int = histogram([f_11_up_5.Ein],binEdges,'Normalization','count');
hold on;
h_ppv = histogram([preEscapeEvents_10([vibr_idx_p prim_idx_p]).Eout],binEdges,...
    'Normalization','count');
h_pes = histogram([preEscapeEvents_10.Eout],binEdges,...
    'Normalization','count');
set(gca,'YScale','Log')
xlabel('KE (eV)');
ylabel('Counts');
legend('Internal','Emission-Optical','Emission');

subplot(1,2,2);
hold off
plot(h_int.Values,h_ppv.Values,'x')
hold on
plot(h_int.Values,h_pes.Values,'o')
set(gca,'XScale','Log')
set(gca,'YScale','Log')
xlabel('Interal Spectrum Count');
ylabel('Counts');
legend('Emissions from optical events removed','All emission events');

%% Investigate if power law is a good description
log_int = log(h_int.Values);
log_ppv = log(h_ppv.Values);
log_pes = log(h_pes.Values);

X = [ones(length(log_int),1) log_int'];

model_pes = fitlm(log_pes,log_int);
model_ppv = fitlm(log_ppv(~isinf(log_ppv)),...
    log_int(~isinf(log_ppv)));
figure(7209);
subplot(1,2,1);
plot(model_pes);
subplot(1,2,2);
plot(model_ppv);

%% Final test
figure(7210);
subplot(1,2,1)
figure(7210);
hold off
semilogy(binCenter,exp(model_pes.Fitted))
hold on
semilogy(binCenter,h_pes.Values)
semilogy(binCenter,h_int.Values)
xlabel('KE(eV)')
ylabel('Counts')
title('Energy Spectra');
legend('Power law model (p = 1.24)','Photoemission','Internal Spectrum')

subplot(1,2,2)
hold off
plot(binCenter,h_int.Values./exp(model_pes.Fitted'));
hold on
plot(binCenter,h_int.Values./h_pes.Values);
legend('Internal/Model_{11}^{1.24}','Internal/Emission')
title('Proportionality');
save('model_pes')

%% Save figures before verification with old data
saveAllFigures('C:\Users\jhansonma\Desktop')

%% Unload the energy scan archive Load the old data
%clear('energyScanArchive');
outputParent    =   strcat('..\\..\\..\\..\\JonathanCodeIO_CXRO\\',...
            'ElectronInteractions\\LEEMRes\\');
outputFolder    =   '20190108_PES_NewVibr';
outputBasePath  =   strcat(outputParent,outputFolder,'\\');

PE_Jan09_Archive = loadEnergyArchive(outputBasePath,'ScanArchive_85eV');

%% The crossing events
f_09_up_5       =   extractZEvent(PE_Jan09_Archive,-5,+1);
f_09_escape     =   extractEscapeEvents(PE_Jan09_Archive);

%% Obtain the histograms
figure(7211);
binEdges = 4:2:88;
subplot(1,2,1);
hold off
h_pes_09 = histogram([f_09_escape.Ein],binEdges);
hold on
h_int_09 = histogram([f_09_up_5.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts');
legend('Photoemission','Internal Spectrum');
title('Jan 09 (Energy dependent MFP) simulation');
binEdges = 4:4:88;

%% 
subplot(1,2,2)
loglog(h_pes_09.Values,h_int_09.Values,'x')
title('Jan 09 (Energy dependent MFP) simulation');
xlabel('Photoemission Counts');
ylabel('Internal Spectrum Counts');

%% Test the power law
pes2int_09 = predict(model_pes,log(h_pes_09.Values)');
figure(7212);
subplot(1,2,1);
hold off
plot(binCenter, h_pes_09.Values);
hold on
plot(binCenter, h_int_09.Values);
plot(binCenter, exp(pes2int_09));
title({'Energy Spectra';...
    'Regression: Jan11 (const IMFP) -> Prediction Jan 09'});
legend('Photoemission','Internal','Model_{11}^{1.24}(Photoemission)')

subplot(1,2,2)
hold off
plot(binCenter, h_int_09.Values./h_pes_09.Values);
hold on
plot(binCenter, h_int_09.Values./exp(pes2int_09)');
legend('Internal/Photoemission','Internal/Model_{11}^{1.24}(Photoemission)')
title('Proportionality')

%% Need a new model
model_pes_09 = fitlm(log(h_pes_09.Values),log(h_int_09.Values));

%% Test the new model
pes2int_09_native = predict(model_pes_09,log(h_pes_09.Values)');
figure(7213);
subplot(1,2,1);
hold off
plot(binCenter, h_pes_09.Values);
hold on
plot(binCenter, h_int_09.Values);
plot(binCenter, exp(pes2int_09_native));
title({'Energy Spectra';...
    'Native model for Jan 09 data'});
legend('Photoemission','Internal','Model_{09}^{1.20}(Photoemission)')
xlabel('KE(eV)')
ylabel('Counts')

subplot(1,2,2)
hold off
plot(binCenter, h_int_09.Values./h_pes_09.Values);
hold on
plot(binCenter, h_int_09.Values./exp(pes2int_09_native)');
plot(binCenter, h_int_09.Values./exp(pes2int_09)');
legend('Internal/Photoemission',...
    'Internal/Model_{09}^{1.20}(Photoemission)',...
    'Internal/Model_{11}^{1.24}(Photoemission)')
title('Proportionality')
xlabel('KE(eV)')
%% The fit for reference
figure(7214)
subplot(1,2,1);
plot(model_pes_09)
title('Log Internal Spectrum (y) vs Log Photoemission count (x1)')

%% Split the data by imfp and see how the data behave
figure(7214);
subplot(1,2,2);
hold off;
loglog(h_pes_09.Values(1:3),h_int_09.Values(1:3),'o')
hold on;
loglog(h_pes_09.Values(4:end),h_int_09.Values(4:end),'x')
xlabel('Photoemission Count');
ylabel('Internal Spectrum Count');
title('Segmentation of data');
legend('Below 16 eV','Above 16 eV','Location','Northwest');
% Doesn't make much sense