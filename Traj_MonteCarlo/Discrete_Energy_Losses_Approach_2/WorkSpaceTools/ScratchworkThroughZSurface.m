%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script investigate the internal electron energy spectra by the means
% of taking statistics of the events that passes through a z contour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Jan 09 archive
addpath('..\..\CommonComponents')
PE_Jan09_Archive    =   loadEnergyArchive('..\\..\\..\\..\\..\\JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\20190108_PES_NewVibr\',...
    'ScanArchive_85eV');
set(0,'DefaultFigureWindowStyle','Docked');
escapeEvents_09 = extractZEvent(PE_Jan09_Archive,0,+1);

%% Jan 09: Event through -5, 

% upwards
[ps5pEV,ps5pIncEV] = extractZEvent(PE_Jan09_Archive,-5,+1);
% downwards
[ps5nEV,ps5nIncEV] = extractZEvent(PE_Jan09_Archive,-5,-1);
% both ways, as a sanity check
[ps5tEV,ps5tIncEV] = extractZEvent(PE_Jan09_Archive,-5,0);
% upward and downwad going electrons should add to to the two-way total
disp(length(ps5pEV)+length(ps5nEV)-length(ps5tEV));

%% Energy spectrum through -5
figure(8000);
binEdges = 3:2:91;
hold off
hx0 = histogram([ps5tEV.Ein],binEdges,'Normalization','count');
hold on;
hpe = histogram([escapeEvents_09.Ein],binEdges,'Normalization','count');
set(gca,'YScale','Log');
xlabel('KE(eV)');
ylabel('Counts')
title({'Simulated electron energy spectra';'17,600 incidences; 1,771,061 events'});
legend('Passage through z = -5 nm (34983 events)','Photoemission (6659 events)');
figure(8001);
binEdges = 4:2:92;
binCneter = (binEdges(1:end-1)+binEdges(2:end))/2;
hold off
hx0 = histogram([ps5tEV.Ein],binEdges,'Normalization','pdf');
hold on;
hpe = histogram([escapeEvents_09.Ein],binEdges,'Normalization','pdf');
set(gca,'YScale','Linear');
xlabel('KE(eV)');
ylabel('Counts')
title({'Simulated electron energy spectra';'17,600 incidences; 1,771,061 events'});
legend('Passage through z = -5 nm (34983 events)','Photoemission (6659 events)');

% Some normalization
figure(8002);
hold off
plot(binCenter,hx0.Values/hx0.Values(41))
hold on;
plot(binCenter,hpe.Values/hpe.Values(41))
plot(binCenter,hpe.Values./(sqrt(binCenter))/(hpe.Values(41)/sqrt(85)))
set(gca,'YScale','Log');
xlabel('KE(eV)');
ylabel('Counts')
title({'Simulated electron energy spectra';'17,600 incidences; 1,771,061 events';...
    'Primary peak normalized'});
legend('Passage through z = -5 nm (34983 events)',...
    'Photoemission (6659 events)',...
    'Photoemission divided by KE^{1/2}');


figure(8010);
hold off
hxp = histogram([ps5pEV.Ein]);
hold on
hxn = histogram([ps5nEV.Ein]);
xlabel('KE(eV)');
ylabel('Counts')
title({'Simulated electron energy spectra';...
    '17,600 incidences; 1,771,061 events';...
    'Passage through z = -5 nm (34983 events)'});
legend('Upward going : 17,680 events','Downward going : 17,303 events');
set(gca,'YScale','Log');

%% Extend it to 1.25, 2.5 7.5, 10 and 10.25 nm
tic
[ps1p25pEV,~] = extractZEvent(PE_Jan09_Archive,-1.25,+1);
[ps1p25nEV,~] = extractZEvent(PE_Jan09_Archive,-1.25,-1);
toc
[ps2p50pEV,~] = extractZEvent(PE_Jan09_Archive,-2.5,+1);
[ps2p50nEV,~] = extractZEvent(PE_Jan09_Archive,-2.5,-1);
toc
[ps7p50tEV,~] = extractZEvent(PE_Jan09_Archive,-7.5,0);
toc
[ps10pEV,~] = extractZEvent(PE_Jan09_Archive,-10,+1);
[ps10nEV,~] = extractZEvent(PE_Jan09_Archive,-10,-1);
toc
[ps10p25pEV,~] = extractZEvent(PE_Jan09_Archive,-12.5,+1);
[ps10p25nEV,~] = extractZEvent(PE_Jan09_Archive,-12.5,-1);
toc

%% Observe the traffic from top to middle
figure(8020);

subplot(2,2,1);
hold off
histogram([ps10pEV.Ein]);
hold on
histogram([ps10nEV.Ein]);
xlabel('KE(eV)');
ylabel('Counts')
title({'Passage through z = -10 nm (17746 events)'});
legend('Upward going','Downward going');
set(gca,'YScale','Log');

subplot(2,2,2);
hold off
histogram([ps5pEV.Ein]);
hold on
histogram([ps5nEV.Ein]);
xlabel('KE(eV)');
ylabel('Counts')
title({'Passage through z = -5 nm (34983 events)'});
set(gca,'YScale','Log');

subplot(2,2,3);
hold off
histogram([ps2p50pEV.Ein]);
hold on
histogram([ps2p50nEV.Ein]);
xlabel('KE(eV)');
ylabel('Counts')
title({'Passage through z = -2.5 nm (30428 events)'});
set(gca,'YScale','Log');


subplot(2,2,4);
hold off
histogram([ps1p25pEV.Ein]);
hold on
histogram([ps1p25nEV.Ein]);
xlabel('KE(eV)');
ylabel('Counts')
title({'Passage through z = -1.25 nm (22602 events)'});
set(gca,'YScale','Log');

%% Let's do the thickness dependence traffic analysis closer to the surface
tic
[ps0p25pEV,~] = extractZEvent(PE_Jan09_Archive,-0.25,+1);
[ps0p25nEV,~] = extractZEvent(PE_Jan09_Archive,-0.25,-1);
toc

tic
[ps0p50pEV,~] = extractZEvent(PE_Jan09_Archive,-0.50,+1);
[ps0p50nEV,~] = extractZEvent(PE_Jan09_Archive,-0.50,-1);
toc

%% Depth dependent traffic round II
figure(8021);
binEdges = 4:2:92;

subplot(2,2,1);
hold off
histogram([ps0p25pEV.Ein],binEdges);
hold on
histogram([ps0p25nEV.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts')
title({'Through z = -0.25 nm';' (10648 events)'});
legend('Upward going','Downward going');
set(gca,'YScale','Linear');

subplot(2,2,2);
hold off
histogram([ps0p50pEV.Ein],binEdges);
hold on
histogram([ps0p50nEV.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts')
title({'Through z = -0.5 nm';' (20363 events)'});
set(gca,'YScale','Linear');

subplot(2,2,3);
hold off
histogram([ps1p25pEV.Ein],binEdges);
hold on
histogram([ps1p25nEV.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts')
title({'Through z = -1.25 nm';' (22602 events)'});
set(gca,'YScale','Linear');

subplot(2,2,4);
hold off
histogram([ps2p50pEV.Ein],binEdges);
hold on
histogram([ps2p50nEV.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts')
title({'Through z = -2.5 nm';' (30428 events)'});
set(gca,'YScale','Linear');

%% Wait, if I put z at zero does that give me the photoemission spectrum
[ps0pEV,~] = extractZEvent(PE_Jan09_Archive,0,+1);
[ps0nEV,~] = extractZEvent(PE_Jan09_Archive,0,-1);
%% Do the numbers add up
disp(length(escapeEvents_09)==length(ps0pEV));
disp(isempty(ps0nEV));
figure(8003);
subplot(2,1,1)
histogram([escapeEvents_09.Ein])
title('Photoemision spectrum')
subplot(2,1,2)
histogram([ps0pEV.Ein])
title('Electrons crossing the resist-vacuum interface')
xlabel('KE (eV)');

%% Try 0.1 nm
[ps0p1pEV,~] = extractZEvent(PE_Jan09_Archive,-0.1,+1);
[ps0p1nEV,~] = extractZEvent(PE_Jan09_Archive,-0.1,-1);

%% Depth dependent traffic round III
figure(8022);
binEdges = 4:2:92;
yscale = 'Linear';

subplot(2,2,1);
hold off
histogram([ps0p1pEV.Ein],binEdges);
hold on
histogram([ps0p1nEV.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts')
title({'Through z = -0.1 nm';' (4691 events)'});
set(gca,'YScale',yscale);

subplot(2,2,2);
hold off
histogram([ps0p25pEV.Ein],binEdges);
hold on
histogram([ps0p25nEV.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts')
title({'Through z = -0.25 nm';' (10648 events)'});
legend('Upward going','Downward going');
set(gca,'YScale',yscale);

subplot(2,2,3);
hold off
histogram([ps0p50pEV.Ein],binEdges);
hold on
histogram([ps0p50nEV.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts')
title({'Through z = -0.5 nm';' (20363 events)'});
set(gca,'YScale',yscale);


subplot(2,2,4);
hold off
histogram([ps2p50pEV.Ein],binEdges);
hold on
histogram([ps2p50nEV.Ein],binEdges);
xlabel('KE(eV)');
ylabel('Counts')
title({'Through z = -2.5 nm';' (30428 events)'});
set(gca,'YScale',yscale)

%% Depth dependent traffic summary


%% Check the angluar distribution before doing directional weighing
figure(8030)
histogram([ps5tEV.theta_in],'Normalization','pdf');
title({'\theta distribution of event crossing z = -5',...
    'Expected functional form : \sim |sin\theta cos\theta|'})
xlabel('\theta in radian');
ylabel('P(\theta)');

%% Division by imfp
% The imfp array
binEdges = 4:2:92;
binCenter = (binEdges(1:end-1)+binEdges(2:end))/2;

imfp_2 = interp1(energyScale,imfp,binCenter);

figure(9000)
hx0 = histogram([ps5tEV.Ein],binEdges,'Normalization','count');
hold on
hpe = histogram([escapeEvents_09.Ein],binEdges,'Normalization','count');

figure(9001);
normSpot = 41;
hold off
plot(binCenter,hx0.Values/hx0.Values(normSpot))
hold on;
plot(binCenter,hpe.Values/hpe.Values(normSpot))
%plot(binCenter,hpe.Values./sqrt(binCenter)/(hpe.Values(normSpot)/sqrt(binCenter(normSpot))))
plot(binCenter,hpe.Values./imfp_2/(hpe.Values(normSpot)/imfp_2(normSpot)))
set(gca,'YScale','Linear');
xlabel('KE(eV)');
ylabel('Normalized Counts')
title({'Simulated electron energy spectra';'17,600 incidences; 1,771,061 events';...
    'Primary peak normalized'});
legend('Passage through z = -5 nm (34983 events)',...
    'Photoemission (6659 events)',...'Photoemission divided by KE^{1/2}',...
    'Photoemission divided by imfp');

%% Finish the depth dependent flux analysis

binEdges = 4:3:92;
binCenter = (binEdges(1:end-1) + binEdges(2:end))/2;

ps0pEV = escapeEvents_09;

figure(666);
hold off
h_up_0 = histogram([ps0pEV.Ein],binEdges,'Normalization','count');
hold on;
h_up_0p10 = histogram([ps0p1pEV.Ein],binEdges,'Normalization','count');
h_up_0p25 = histogram([ps0p25pEV.Ein],binEdges,'Normalization','count');
h_up_0p50 = histogram([ps0p50pEV.Ein],binEdges,'Normalization','count');
h_up_2p50 = histogram([ps2p50pEV.Ein],binEdges,'Normalization','count');

figure(8023);
marker = '-';
[~,normp] = min(abs(binCenter-85));
hold off
plot(binCenter,h_up_0.Values./h_up_0.Values(normp+1),marker)
hold on
plot(binCenter,h_up_0p10.Values./h_up_0p10.Values(normp+1),marker)
plot(binCenter,h_up_0p25.Values./h_up_0p25.Values(normp+1),marker)
plot(binCenter,h_up_0p50.Values./h_up_0p50.Values(normp+1),marker)
plot(binCenter,h_up_2p50.Values./h_up_2p50.Values(normp+1),marker)
set(gca,'YScale','Linear');
legend({'Photoemission';...
    '-0.1 nm';'-0.25 nm';'-0.5 nm';'-2.5 nm'});
title({'Simulated upward flux at different depths';...
    'Energy Probability Density'});
xlabel('KE (eV)')
ylabel('Normalized Counts')
set(gca,'YScale','Linear')

% Downward flux
figure(666);
hold off
h_dn_0p10 = histogram([ps0p1nEV.Ein],binEdges,'Normalization','count');
hold on;
h_dn_0p25 = histogram([ps0p25nEV.Ein],binEdges,'Normalization','count');
h_dn_0p50 = histogram([ps0p50nEV.Ein],binEdges,'Normalization','count');
h_dn_2p50 = histogram([ps2p50nEV.Ein],binEdges,'Normalization','count');

figure(8024);
[~,normp] = min(abs(binCenter-85));
marker = '-';
hold off
plot(binCenter,h_dn_0p10.Values./h_dn_0p10.Values(normp+1),marker)
hold on
plot(binCenter,h_dn_0p25.Values./h_dn_0p25.Values(normp+1),marker)
plot(binCenter,h_dn_0p50.Values./h_dn_0p50.Values(normp+1),marker)
plot(binCenter,h_dn_2p50.Values./h_dn_2p50.Values(normp+1),marker)
set(gca,'YScale','Linear');
legend({'-0.1 nm';'-0.25 nm';'-0.5 nm';'-2.5 nm'});
title({'Simulated downward flux at different depths';...
    'Energy Probability Density'});
xlabel('KE (eV)');



%% Exporting the graphs
addpath('..\..\CommonComponents');
saveAllFigures(outputBasePath);
