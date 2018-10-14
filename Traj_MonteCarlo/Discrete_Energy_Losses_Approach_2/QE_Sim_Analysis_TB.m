%% Analysis of the QE results
clear;
clc;
close all;

addpath('F:\Documents and Settings\sbhattarai\My Documents\Research\Globals\xticklabel_rotate');

filepath='QE_Stats_low10AllPAGs\';
% filepath='QE_Stats_Sims_6\';
% filepath='QE_Stats_lowEthr\';

% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=91.00eV_pagEa=3_Scatt-Elim=20_PAG-EaMin=5_1024Trials.mat';
% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=30.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_512Trials.mat';
filename='QEStats_Low10AllPAG_Dose=1.00_RcnRad=2.00_PAG=0.40_E=91.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_1024Trials.mat';
filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=91.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_1024Trials_ReactThenMove_BKUP.mat';
filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=91.00eV_pagEa=5_Scatt-Elim=20_PAG-EaMin=5_1024Trials_RcnThenMove.mat';
% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=91.00eV_pagEa=5_Scatt-Elim=18.00_PAG-EaMin=5_3500Trials_RcnMove.mat';
data91=load([filepath filename]);
% data30=load([filepath filename]);

% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=80.00eV_pagEa=3_Scatt-Elim=20_PAG-EaMin=5_1024Trials.mat';
% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=40.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_512Trials.mat';
filename='QEStats_Low10AllPAG_Dose=1.00_RcnRad=2.00_PAG=0.40_E=80.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_1024Trials.mat';
filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=80.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_1024Trials_ReactThenMove_BKUP.mat';
filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=80.00eV_pagEa=5_Scatt-Elim=20_PAG-EaMin=5_1024Trials_RcnThenMove.mat';
% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=80.00eV_pagEa=5_Scatt-Elim=18.00_PAG-EaMin=5_3500Trials_RcnMove.mat';
data80=load([filepath filename]);
% data40=load([filepath filename]);

filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=50.00eV_pagEa=3_Scatt-Elim=20_PAG-EaMin=5_1024Trials.mat';
% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=50.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_512Trials.mat';
filename='QEStats_Low10AllPAG_Dose=1.00_RcnRad=2.00_PAG=0.40_E=50.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_1024Trials.mat';
filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=50.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_1024Trials_ReactThenMove.mat';
filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=50.00eV_pagEa=5_Scatt-Elim=20_PAG-EaMin=5_1024Trials_RcnThenMove.mat';
% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=50.00eV_pagEa=5_Scatt-Elim=18.00_PAG-EaMin=5_3500Trials_RcnMove.mat';
data50=load([filepath filename]);
% data50=load([filepath filename]);

filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=30.00eV_pagEa=3_Scatt-Elim=20_PAG-EaMin=5_1024Trials.mat';
% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=60.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_512Trials.mat';
filename='QEStats_Low10AllPAG_Dose=1.00_RcnRad=2.00_PAG=0.40_E=30.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_1024Trials.mat';
filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=30.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_1024Trials_ReactThenMove.mat';
filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=30.00eV_pagEa=5_Scatt-Elim=20_PAG-EaMin=5_1024Trials_RcnThenMove.mat';
% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=30.00eV_pagEa=5_Scatt-Elim=18.00_PAG-EaMin=5_3500Trials_RcnMove.mat';
data30=load([filepath filename]);
% data60=load([filepath filename]);

% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=80.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_512Trials.mat';
% filename='QEStats_Dose=1.00_RcnRad=2.00_PAG=0.40_E=80.00eV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_512Trials.mat';
% data80=load([filepath filename]);

%% Add up all of polymers and PAG excitations into 1 matrix [wrong, need the 6eV-polym matrix which I didn't save]
clc;
close all;

data30_2={};
data50_2={};
data80_2={};
data91_2={};

data30_2.acidimg_global=data30.acidimg_global+data30.ionimg_global;
idx=find(data30_2.acidimg_global>0);
data30_2.xyz_acids=[data30.univ.grid.x(idx) data30.univ.grid.y(idx) data30.univ.grid.z(idx)];
data30_2.radius_acids_global=sqrt((data30_2.xyz_acids(:,1)-data30.xval).^2+(data30_2.xyz_acids(:,2)-data30.yval).^2+(data30_2.xyz_acids(:,3)-data30.zval).^2);
data30_2.nacids=sum(data30_2.acidimg_global(:))./data30.ntrials;
data30_2.meanacids=mean([data30.nacids_thrutrial+data30.nions_thrutrial]);
data30_2.stdacids=std([data30.nacids_thrutrial+data30.nions_thrutrial]);

data50_2.acidimg_global=data50.acidimg_global+data50.ionimg_global;
idx=find(data50_2.acidimg_global>0);
data50_2.xyz_acids=[data50.univ.grid.x(idx) data50.univ.grid.y(idx) data50.univ.grid.z(idx)];
data50_2.radius_acids_global=sqrt((data50_2.xyz_acids(:,1)-data50.xval).^2+(data50_2.xyz_acids(:,2)-data50.yval).^2+(data50_2.xyz_acids(:,3)-data50.zval).^2);
data50_2.nacids=sum(data50_2.acidimg_global(:))./data50.ntrials;
data50_2.meanacids=mean([data50.nacids_thrutrial+data50.nions_thrutrial]);
data50_2.stdacids=std([data50.nacids_thrutrial+data50.nions_thrutrial]);

data80_2.acidimg_global=data80.acidimg_global+data80.ionimg_global;
idx=find(data80_2.acidimg_global>0);
data80_2.xyz_acids=[data80.univ.grid.x(idx) data80.univ.grid.y(idx) data80.univ.grid.z(idx)];
data80_2.radius_acids_global=sqrt((data80_2.xyz_acids(:,1)-data80.xval).^2+(data80_2.xyz_acids(:,2)-data80.yval).^2+(data80_2.xyz_acids(:,3)-data80.zval).^2);
data80_2.nacids=sum(data80_2.acidimg_global(:))./data80.ntrials;
data80_2.meanacids=mean([data80.nacids_thrutrial+data80.nions_thrutrial]);
data80_2.stdacids=std([data80.nacids_thrutrial+data80.nions_thrutrial]);

data91_2.acidimg_global=data91.acidimg_global+data91.ionimg_global;
idx=find(data91_2.acidimg_global>0);
data91_2.xyz_acids=[data91.univ.grid.x(idx) data91.univ.grid.y(idx) data91.univ.grid.z(idx)];
data91_2.radius_acids_global=sqrt((xyz_acids(:,1)-data91.xval).^2+(xyz_acids(:,2)-data91.yval).^2+(xyz_acids(:,3)-data91.zval).^2);
data91_2.nacids=sum(data91_2.acidimg_global(:))./data91.ntrials;
data91_2.meanacids=mean([data91.nacids_thrutrial+data91.nions_thrutrial]);
data91_2.stdacids=std([data91.nacids_thrutrial+data91.nions_thrutrial]);

figure;
plot([30 50 80 91],[data30_2.nacids data50_2.nacids data80_2.nacids data91_2.nacids],'-ob');

figure;
plot([30 50 80 91],[data30_2.nacids data50_2.nacids data80_2.nacids data91_2.nacids]./([30 50 80 91]),'-ob');

figure;
errorbar([30 50 80 91],[data30_2.meanacids data50_2.meanacids data80_2.meanacids data91_2.meanacids],[data30_2.stdacids data50_2.stdacids data80_2.stdacids data91_2.stdacids],'-ob');

% figure;
% errorbar([30 50 80 91],[data30_2.meanacids data50_2.meanacids data80_2.meanacids data91_2.meanacids],[data30_2.stdacids data50_2.stdacids data80_2.stdacids data91_2.stdacids],'-ob');

%% Get mean and stdevs of the acids
clc
close all;

gcf_pos=[106   282   845   642];
gca_pos=[ 0.1773    0.1310    0.7277    0.7940];

E=[30 50 80 91];
% meanacids=[mean(data30.nacids_thrutrial) mean(data40.nacids_thrutrial) mean(data50.nacids_thrutrial) mean(data60.nacids_thrutrial) mean(data80.nacids_thrutrial)];
meanacids=[mean(data30.nacids_thrutrial) mean(data50.nacids_thrutrial) mean(data80.nacids_thrutrial) mean(data91.nacids_thrutrial)];
% meanacids=[sum(data30.nacids_thrutrial) sum(data50.nacids_thrutrial) sum(data80.nacids_thrutrial) sum(data91.nacids_thrutrial)];

% stdacids=[std(data30.nacids_thrutrial) std(data40.nacids_thrutrial) std(data50.nacids_thrutrial) std(data60.nacids_thrutrial) std(data80.nacids_thrutrial)];
stdacids=[std(data30.nacids_thrutrial) std(data50.nacids_thrutrial) std(data80.nacids_thrutrial) std(data91.nacids_thrutrial)];

[n30,x30]=hist(data30.nacids_thrutrial,100);    %x30=x30(n30~=0);    n30=n30(n30~=0);    n30=n30./sum(n30);
n30=n30./sum(n30);
% [n40,x40]=hist(data40.nacids_thrutrial,100);    x40=x40(n40~=0);    n40=n40(n40~=0);    n40=n40./sum(n40);
[n50,x50]=hist(data50.nacids_thrutrial,100);    %x50=x50(n50~=0);    n50=n50(n50~=0);    n50=n50./sum(n50);
n50=n50./sum(n50);
% [n60,x60]=hist(data60.nacids_thrutrial,100);    x60=x60(n60~=0);    n60=n60(n60~=0);    n60=n60./sum(n60);
[n80,x80]=hist(data80.nacids_thrutrial,100);    %x80=x80(n80~=0);    n80=n80(n80~=0);    n80=n80./sum(n80);
n80=n80./sum(n80);
[n91,x91]=hist(data91.nacids_thrutrial,100);    %x91=x91(n91~=0);    n91=n91(n91~=0);    n91=n91./sum(n91);
n91=n91./sum(n91);

figure;
plot(E,meanacids,'-ob','linewidth',3.0,'MarkerSIze',12);
% errorbar(E,meanacids,stdacids,'-ob','linewidth',3.0,'MarkerSIze',12);
xlabel('Energy (eV)','fontsize',30);
ylabel('# of Acids','fontsize',30);
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
% ylim([0 2]);

figure;
plot(E,meanacids./meanacids(end),'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
plot(E,E./E(end),'--b','linewidth',3.0);
xlabel('Energy (eV)','fontsize',30);
ylabel('# of Acids/Energy','fontsize',30);
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);


figure;
plot(E,meanacids./E,'-ob','linewidth',3.0,'MarkerSize',12);
% hold on;
% plot(E,E./E(end),'--b','linewidth',3.0);
xlabel('Energy (eV)','fontsize',30);
ylabel('Acids per eV','fontsize',30);
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
gcf_pos=[548   562   499   420];
gca_pos=[0.1300    0.1100    0.7750    0.7329];

figure;
bar(x30,n30,1);
stem(x30(n30~=0),n30(n30~=0),'-ob','linewidth',3.0,'MarkerSize',12);
set(gca,'fontsize',30,'linewidth',3.0);box on;
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
title('30 eV');
xlim([-0.1 8]);ylim([0 0.6]);

figure;
bar(x50,n50,1);
stem(x50(n50~=0),n50(n50~=0),'-ob','linewidth',3.0,'MarkerSize',12);
set(gca,'fontsize',30,'linewidth',3.0);box on;
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
title('50 eV');
xlim([-0.1 8]);ylim([0 0.6]);

figure;
bar(x80,n80,1);
stem(x80(n80~=0),n80(n80~=0),'-ob','linewidth',3.0,'MarkerSize',12);
set(gca,'fontsize',30,'linewidth',3.0);box on;
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
title('80 eV');
xlim([-0.1 8]);ylim([0 0.6]);

figure;
bar(x91,n91,1);
stem(x91(n91~=0),n91(n91~=0),'-ob','linewidth',3.0,'MarkerSize',12);
set(gca,'fontsize',30,'linewidth',3.0);box on;
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
title('91 eV');
xlim([-0.1 8]);ylim([0 0.6]);

%% Get the activation ID
clc;
% close all;

%%% results from running the cell previously:
Etmp=[30;50;80;91];
Ntot=[608 104 451; % Ntot without the wasted energy
      970 179 1103;
      1398 313 2149;
      1524 387 2497]; % [6eV-polym Acid SE]

Ntot=[608 574 451 104;
      970 942 1103 179;
      1398 1253 2149 313;
      1524 1481 2497 387]; % [6eV-Polym Eloss<pagEaMin SE Acid]

Ntot=[863 812 193 133;
      912 926 1054 186;
      1367 1328 2071 284;
      1526 1516 2431 368]; % [6eV-Polym Eloss<pagEaMin SE Acid]

Ntot_perelec=Ntot./1024;
acids_perelec=sum(Ntot_perelec(:,1:2),2);

Ntot_perelec_perE=Ntot_perelec./repmat(Etmp,1,size(Ntot,2));
acids_perelec_perE=acids_perelec./Etmp;

gcf_pos2=[534   207   921   721];
gca_pos2=[ 0.1596    0.1137    0.7660    0.8419];

figure;
plot(Etmp(:,1),Ntot_perelec(:,1),'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
% plot(Etmp(:,1),Ntot_perelec(:,2),'-ok','linewidth',3.0,'MarkerSize',12);
% plot(Etmp(:,1),Ntot_perelec(:,3),'-or','linewidth',3.0,'MarkerSize',12);
plot(Etmp(:,1),Ntot_perelec(:,4),'-ok','linewidth',3.0,'MarkerSize',12);
plot(Etmp(:,1),Ntot_perelec(:,1)+Ntot_perelec(:,4),'-or','linewidth',3.0,'MarkerSize',12);
xlabel('Energy (eV)','fontsize',15);
ylabel('Events/Electron','fontsize',15);
set(gca,'fontsize',30,'linewidth',3.0);box on;
set(gca,'Position',gca_pos2);
set(gcf,'Position',gcf_pos2);
legend('6eV-Polym','Eloss<5eV','SE','Acid');
legend('Polym','Acid','Polym+Acid');


figure;
plot(Etmp(:,1),acids_perelec,'-ob','linewidth',3.0,'MarkerSize',12);
xlabel('Energy (eV)','fontsize',15);
ylabel('Acids/Electron','fontsize',15);
set(gca,'fontsize',30,'linewidth',3.0);box on;
set(gca,'Position',gca_pos2);
set(gcf,'Position',gcf_pos2);
% legend('6eV-Polym','Eloss<5eV','SE','Acid');

figure;
plot(Etmp(:,1),acids_perelec_perE,'-ob','linewidth',3.0,'MarkerSize',12);
xlabel('Energy (eV)','fontsize',15);
ylabel('Acids/Electron/eV','fontsize',15);
set(gca,'fontsize',30,'linewidth',3.0);box on;
set(gca,'Position',gca_pos2);
set(gcf,'Position',gcf_pos2);
% legend('6eV-Polym','Eloss<5eV','SE','Acid');

%%
figure;
plot(Etmp(:,1),Ntot_perelec_perE(:,1),'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
plot(Etmp(:,1),Ntot_perelec_perE(:,2),'-ok','linewidth',3.0,'MarkerSize',12);
plot(Etmp(:,1),Ntot_perelec_perE(:,3),'-or','linewidth',3.0,'MarkerSize',12);
plot(Etmp(:,1),Ntot_perelec_perE(:,4),'-og','linewidth',3.0,'MarkerSize',12);
xlabel('Energy (eV)','fontsize',15);
ylabel('Events/Electron/eV','fontsize',15);
set(gca,'fontsize',30,'linewidth',3.0);box on;
set(gca,'Position',gca_pos2);
set(gcf,'Position',gcf_pos2);
legend('6eV-Polym','Eloss<5eV','SE','Acid');

%%% Acid fraction
figure;
plot(Etmp(:,1),Ntot_perelec_perE(:,4)./sum(Ntot_perelec_perE,2),'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
% plot(Etmp(:,1),Ntot_perelec_perE(:,2),'-ok','linewidth',3.0,'MarkerSize',12);
% plot(Etmp(:,1),Ntot_perelec_perE(:,3),'-or','linewidth',3.0,'MarkerSize',12);
xlabel('Energy (eV)','fontsize',15);
ylabel('Acid/Polym','fontsize',15);
set(gca,'fontsize',30,'linewidth',3.0);box on;
set(gca,'Position',gca_pos2);
set(gcf,'Position',gcf_pos2);
% legend('6eV-Polym','Acid','SE');
ylim([0 0.1]);

%%% SE fraction
figure;
plot(Etmp(:,1),Ntot_perelec(:,3)./sum(Ntot_perelec,2),'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
% plot(Etmp(:,1),Ntot_perelec_perE(:,2),'-ok','linewidth',3.0,'MarkerSize',12);
% plot(Etmp(:,1),Ntot_perelec_perE(:,3),'-or','linewidth',3.0,'MarkerSize',12);
xlabel('Energy (eV)','fontsize',15);
ylabel('SE/Total','fontsize',15);
set(gca,'fontsize',30,'linewidth',3.0);box on;
set(gca,'Position',gca_pos2);
set(gcf,'Position',gcf_pos2);
% legend('6eV-Polym','Acid','SE');
% ylim([0 0.1]);

%% Just read the activation data
close all;
clc;

gcf_pos=[680         161        1191         817];
gca_pos=[0.1643    0.3133    0.7266    0.6117];
act_global=data91.act_global;

index = find(strcmp('Eloss<pagEaMin', act_global), Inf);
% act_global(index)=[];

[act_global_unique,idx1,idx2]=unique(act_global);

N=[];
act_tmp={};
for i = 1:length(idx1)
    idx=find(idx2==i);
    act_tmp{i}=act_global(idx);
    N(i)=length(idx);
end
N
act_global_unique

fprintf('sum(N) = %.2f\n',sum(N));
fprintf('sum(N/1024) = %.2f\n',sum(N./1024));
fprintf('****************\n');

figure;
bar(N);
ax=gca;
% set(ax,'XTickLabel',{act_global_unique{1}, act_global_unique{2}, act_global_unique{3}});
set(ax,'XTickLabel',{act_global_unique{1}, act_global_unique{2}, act_global_unique{3}, act_global_unique{4}});
set(ax,'fontsize',30);
set(ax,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
xticklabel_rotate([],45,[],'fontsize',30);
set(gca,'linewidth',3.0);


figure;
bar(N./1024);
ax=gca;
% set(ax,'XTickLabel',{act_global_unique{1}, act_global_unique{2}, act_global_unique{3}});
set(ax,'XTickLabel',{act_global_unique{1}, act_global_unique{2}, act_global_unique{3}, act_global_unique{4}});
set(ax,'fontsize',30);
set(ax,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
xticklabel_rotate([],45,[],'fontsize',30);
set(gca,'linewidth',3.0);

%% Compare ELF with analytical:
close all;

optdata_path='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';
% scattdata.optical=load([optdata_path 'Sp_IMFP_Inelastic_Components_Ef=0p5eV_Elossmin=0.001eV_Erange=[5,1000]_DDCSData.mat']); % _v2 has data structures better suited to the subsequently called programs
% scattdata.optical=load([optdata_path 'Sp_withICSData_IMFP_Inelastic_Components_Ef=10eV_Elossmin=0.001eV_Erange=[5,1000]_DDCSData.mat']);
% scattdata.optical=load([optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData_EQCons=Pines.mat']);
scattdata.optical=load([optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData.mat']);

pathname='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';
filename='DDCSdata_Ef=0p5_Elossmin=0.001eV_Erange=[5,1000].mat';
filename='DDCSdata_withICSData_Ef=10_Elossmin=0.001eV_Erange=[5,1000].mat';
filename='DDCSdata_Fuji_Ef=15.5_Elossmin=3eV_Erange=[19,200]_EQCons=Pines.mat';
filename='DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,200].mat';
scattdata.optical.inel_dcsdata=load([pathname filename]);

% [30,50,80,92]: [8,18,33,39]

xvec=scattdata.optical.inel_dcsdata.Elossmat(8,:);
yvec=scattdata.optical.inel_dcsdata.dsigdE(8,:);

gcf_pos=[680   256   945   722];
gca_pos=[0.1672    0.1219    0.7378    0.8031];

[n1,x1]=hist(data30.Eloss(:,1),15);
n1=n1./trapz(x1,n1);

[n2,x2]=hist(data50.Eloss(:,1),15);
n2=n2./trapz(x2,n2);

[n3,x3]=hist(data80.Eloss(:,1),15);
n3=n3./trapz(x3,n3);

[n4,x4]=hist(data91.Eloss(:,1),15);
n4=n4./trapz(x4,n4);

% 30eV: idx = 8;
figure;
bar(x1,n1,1);
hold on;
Evec=scattdata.optical.inel_dcsdata.Elossmat(8,:);
pdfvec=scattdata.optical.inel_dcsdata.dsigdE(8,:);
pdfvec=pdfvec./trapz(Evec,pdfvec);
plot(Evec,pdfvec,'r','linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);

% 50 eV: idx = 18;
figure;
bar(x2,n2,1);
hold on;
Evec=scattdata.optical.inel_dcsdata.Elossmat(18,:);
pdfvec=scattdata.optical.inel_dcsdata.dsigdE(18,:);
pdfvec=pdfvec./trapz(Evec,pdfvec);
plot(Evec,pdfvec,'r','linewidth',4.0);
xlabel('Energy loss (eV)','fontsize',30);
ylabel('Prob. Density. Func','fontsize',30);
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);

% 80 eV: idx = 33;
figure;
bar(x3,n3,1);
hold on;
Evec=scattdata.optical.inel_dcsdata.Elossmat(33,:);
pdfvec=scattdata.optical.inel_dcsdata.dsigdE(33,:);
pdfvec=pdfvec./trapz(Evec,pdfvec);
plot(Evec,pdfvec,'r','linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);

% 90 eV: idx = 38;
figure;
bar(x4,n4,1);
hold on;
Evec=scattdata.optical.inel_dcsdata.Elossmat(38,:);
pdfvec=scattdata.optical.inel_dcsdata.dsigdE(38,:);
pdfvec=pdfvec./trapz(Evec,pdfvec);
plot(Evec,pdfvec,'-r','linewidth',6.0);
xlabel('Energy loss (eV)','fontsize',30);
ylabel('Prob. Density. Func','fontsize',30);
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);

%% 3-D distributions that ANdy wanted
clc;
% close all;

analysis_entity=data30;

idx=find(analysis_entity.acidimg_global>0);
[xidx,yidx,zidx]=ind2sub(size(analysis_entity.acidimg),idx);
acidvec=[];
for i = 1:length(xidx)
    acidvec(i)=analysis_entity.acidimg_global(xidx(i),yidx(i),zidx(i));
end

xo=analysis_entity.xval;
yo=analysis_entity.yval;
zo=analysis_entity.zval;

idx1=find(analysis_entity.univ.grid.x==xo & analysis_entity.univ.grid.y==yo & analysis_entity.univ.grid.z==zo);
[xidx1,yidx1,zidx1]=ind2sub(size(analysis_entity.acidimg),idx1);

xvec=xo-5:xo+5;
yvec=yo-5:yo+5;
zvec=zo-5:zo+5;

[xmat,zmat,ymat]=meshgrid(xvec);

r=sqrt(xmat.^2+ymat.^2+zmat.^2);

count=1;
rho=[];
rcnrad=1;
for i = 1:length(xvec)
    for j = 1:length(yvec)
        for k = 1:length(zvec)
            fprintf('Co-ordinate %d of %d\n',count,length(xvec)*length(yvec)*length(zvec));
            [yesno,rho(k,i,j),pagidx,npixels_vol]=pag_avail([xvec(i) yvec(j) zvec(k)],analysis_entity.univ.grid,analysis_entity.acidimg_global,rcnrad);
            rho(k,i,j)=rho(k,i,j)./npixels_vol;
            count=count+1;
        end
    end
end
% rho=rho./(4/3*pi*rcnrad^3);

% figure;
% hist(sqrt((xidx-xidx1).^2+(yidx-yidx1).^2+(zidx-zidx1).^2))

% figure;
% plot(acidvec);

%% density vs. radius: average the density
clc;
close all;

radius=0:10;

[xmat,zmat,ymat]=meshgrid(xvec);
tmpgrid.x=xmat;
tmpgrid.y=ymat;
tmpgrid.z=zmat;

navail2=[];
for i = 1:length(radius)
    if i>=2
        radius_old=radius(i-1);
    else
        radius_old=0;
    end
    [yesno,navail2(i),navail2_std(i),pagidx]=NumatRad([analysis_entity.xval analysis_entity.yval analysis_entity.zval],tmpgrid,rho,[radius_old radius(i)]);
end
% radius=[0 radius];
% navail2=[0 navail2];

% dndr=dydx(radius,navail2);
dndr=diff(navail2)./diff(radius);
dndr=[0 dndr];

figure;
plot(radius,navail2,'-o');

figure;
plot(radius,dndr,'-o');

%% acid gen probability per spherical shell radius
clc;
close all;

radius=0:1:20;

% xo=analysis_entity.xval;
% yo=analysis_entity.yval;
% zo=analysis_entity.zval;
% 
% idx1=find(analysis_entity.univ.grid.x==xo & analysis_entity.univ.grid.y==yo & analysis_entity.univ.grid.z==zo);
% [xidx1,yidx1,zidx1]=ind2sub(size(analysis_entity.acidimg),idx1);

analysis_entity={};
analysis_entity{length(analysis_entity)+1}=data30;
% analysis_entity{length(analysis_entity)+1}=data40;
analysis_entity{length(analysis_entity)+1}=data50;
% analysis_entity{length(analysis_entity)+1}=data60;
analysis_entity{length(analysis_entity)+1}=data80;
analysis_entity{length(analysis_entity)+1}=data91;

navail2=[];
dndr=[];
radius2=[];
navail2=[];
dndr=[];
dndr_norm=[];
imfp_mean=[];
imfp_std=[];
for j = 1:length(analysis_entity)
    navail=[];
    for i = 1:length(radius)
        [yesno,navail(1,i),pagidx,npixels_vol]=pag_avail([analysis_entity{j}.xval analysis_entity{j}.yval analysis_entity{j}.zval],analysis_entity{j}.univ.grid,analysis_entity{j}.acidimg_global,radius(i));
    %     dndr=diff(navail2)./diff(radius);
%         navail(1,i)/npixels_vol;
        
    end
%     navail=navail./1024;
%     navail=navail./analysis_entity{j}.event{1}.Ese;
    
    radius2(j,:)=[radius];
    navail2(j,:)=[navail];
%     dndr_tmp=dydx(radius,navail);
    dndr(j,:)=[0 diff(navail)./diff(radius)];
%     dndr(j,:)=[0 dndr_tmp dndr_tmp(end)];
    dndr_norm(j,:)=dndr(j,:)./max(dndr(j,:));

    xtmp=radius2(j,:);
    pdf_tmp=dndr(j,:);
    pdf_tmp=pdf_tmp./trapz(xtmp,pdf_tmp);
    imfp_mean(j)=trapz(xtmp,xtmp.*pdf_tmp);
    imfp_std(j)=sqrt(trapz(xtmp,xtmp.^2.*pdf_tmp)-(imfp_mean(j))^2);

end

%%% fit an analytical function
meanval=2;
sigma=1;
radius_fit=0:0.1:20;
psf1=1/(sigma*sqrt(2*pi)).*exp(-(radius_fit-meanval).^2./(2*sigma^2));                  psf1=psf1./max(psf1).*max(dndr(1,:));
psf2=1/(sigma*sqrt(2*pi)).*radius_fit.*exp(-(radius_fit-meanval).^2./(2*sigma^2));      psf2=psf2./max(psf2).*max(dndr(1,:));

gcf_pos=[ 268   345   750   559];
gca_pos=[ 0.1300    0.1100    0.7750    0.8150];

figure;
% plot(radius2',navail2','-o','linewidth',3.0);
plot(radius2(1,:),navail2(1,:),'--ob','linewidth',3.0);
hold on;
plot(radius2(2,:),navail2(2,:),'--ok','linewidth',3.0);
plot(radius2(3,:),navail2(3,:),'--or','linewidth',3.0);
plot(radius2(4,:),navail2(4,:),'--og','linewidth',3.0);
xlabel('Radius (nm)');ylabel('# Available');
set(gca,'fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);
set(gca,'Position',gca_pos);

figure;
% plot(radius2',dndr','-o','linewidth',3.0);
plot(radius2(1,:),dndr(1,:),'-ob','linewidth',3.0);
hold on;
plot(radius2(2,:),dndr(2,:),'-ok','linewidth',3.0);
plot(radius2(3,:),dndr(3,:),'-or','linewidth',3.0);
plot(radius2(4,:),dndr(4,:),'-og','linewidth',3.0);
% plot(radius_fit,psf1,'--k','linewidth',3.0);
% plot(radius_fit,psf2,'--k','linewidth',3.0);
set(gca,'fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);
set(gca,'Position',gca_pos);


figure;
plot(radius2',dndr_norm','-o','linewidth',3.0);
plot(radius2(1,:),dndr_norm(1,:),'-ob','linewidth',3.0);
hold on;
plot(radius2(2,:),dndr_norm(2,:),'-ok','linewidth',3.0);
plot(radius2(3,:),dndr_norm(3,:),'-or','linewidth',3.0);
plot(radius2(4,:),dndr_norm(4,:),'-og','linewidth',3.0);
set(gca,'fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);
set(gca,'Position',gca_pos);

%% plot mean and stdev of acids
close all;

gcf_pos=[803   325   804   630];
gca_pos=[0.1418    0.1397    0.7632    0.7853];
figure;
errorbar(E,meanacids,stdacids);
% plot(E,meanacids,'-ob','linewidth',3.0,'MarkerSIze',12);
% plot(E,meanacids./meanacids(end),'-ob','linewidth',3.0,'MarkerSIze',12);
% hold on;plot(E,E./E(end),'--b','linewidth',3.0);
xlabel('Electron Energy (eV)','fontsize',30);
ylabel('Number of Acids','fontsize',30);
set(gca,'fontsize',30,'Position',gca_pos,'linewidth',3.0);
set(gcf,'Position',gcf_pos);

figure;
errorbar(E,meanacids,stdacids);
plot(E,meanacids,'-ob','linewidth',3.0,'MarkerSIze',12);
plot(E,meanacids./meanacids(end),'-ob','linewidth',3.0,'MarkerSIze',12);
% hold on;plot(E,E./E(end),'--b','linewidth',3.0);
xlabel('Electron Energy (eV)','fontsize',30);
ylabel('Number of Acids','fontsize',30);
set(gca,'fontsize',30,'Position',gca_pos,'linewidth',3.0);
set(gcf,'Position',gcf_pos);

figure;
plot(E,meanacids./E,'-ob','linewidth',3.0,'MarkerSize',12);
xlabel('Electron Energy (eV)','fontsize',30);
ylabel('Number of Acids','fontsize',30);
set(gca,'fontsize',30,'Position',gca_pos,'linewidth',3.0);
set(gcf,'Position',gcf_pos);

figure;
plot(E,stdacids,'-ob','linewidth',3.0);
xlabel('Electron Energy (eV)','fontsize',30);
ylabel('Number of Acids','fontsize',30);

%% 2-D spatial distribution views

close all;

gcf_pos=[680   333   810   645];
gca_pos=[0.1300    0.1100    0.6894    0.8150];

refidx=[26 25]; % [z,x] in pixels
% refidx=[26 25 25]; % [row,col,third-dim] in pixels

data30_2=[];
data50_2=[];
data80_2=[];
data91_2=[];

zx=sum(data30.acidimg_global,3);
% zx=data30.acidimg_global;
zx2=zx(zx>0);acids_per_px=sort(zx2);
[x1,n1]=hist(acids_per_px,20);
acids_per_px=unique(acids_per_px);
radvec=[];
Ntmp=[];
for i = 1:length(acids_per_px)
    idx=find(zx==acids_per_px(i));
    [xidx,yidx]=ind2sub(size(zx),idx);
%     [xidx,yidx,zidx]=ind2sub(size(zx),idx);
    rad_tmp=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2);
%     rad_tmp=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2 + (zidx-refidx(3)).^2);
    radvec=[radvec;repmat(rad_tmp,acids_per_px(i),1)];
end

idx=find(zx>0);
[xidx,yidx]=ind2sub(size(zx),idx);
% data30_2.radius_acids_global=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2);
data30_2.radius_acids_global=radvec;

zx=sum(data50.acidimg_global,3);
% zx=data50.acidimg_global;
zx2=zx(zx>0);acids_per_px=sort(zx2);
[x1,n1]=hist(acids_per_px,20);
acids_per_px=unique(acids_per_px);
radvec=[];
for i = 1:length(acids_per_px)
    idx=find(zx==acids_per_px(i));
    [xidx,yidx]=ind2sub(size(zx),idx);
%     [xidx,yidx,zidx]=ind2sub(size(zx),idx);
    rad_tmp=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2);
%     rad_tmp=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2 + (zidx-refidx(3)).^2);
    radvec=[radvec;repmat(rad_tmp,acids_per_px(i),1)];
end
idx=find(zx>0);
[xidx,yidx]=ind2sub(size(zx),idx);
% data50_2.radius_acids_global=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2);
data50_2.radius_acids_global=radvec;

zx=sum(data80.acidimg_global,3);
% zx=data80.acidimg_global;
zx2=zx(zx>0);acids_per_px=sort(zx2);
[x1,n1]=hist(acids_per_px,20);
acids_per_px=unique(acids_per_px);
radvec=[];
for i = 1:length(acids_per_px)
    idx=find(zx==acids_per_px(i));
    [xidx,yidx]=ind2sub(size(zx),idx);
%     [xidx,yidx,zidx]=ind2sub(size(zx),idx);
    rad_tmp=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2);
%     rad_tmp=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2 + (zidx-refidx(3)).^2);
    radvec=[radvec;repmat(rad_tmp,acids_per_px(i),1)];
end
idx=find(zx>0);
[xidx,yidx]=ind2sub(size(zx),idx);
% data80_2.radius_acids_global=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2);
data80_2.radius_acids_global=radvec;

zx=sum(data91.acidimg_global,3);
% zx=data91.acidimg_global;
zx2=zx(zx>0);acids_per_px=sort(zx2);
[x1,n1]=hist(acids_per_px,20);
acids_per_px=unique(acids_per_px);
radvec=[];
for i = 1:length(acids_per_px)
    idx=find(zx==acids_per_px(i));
    [xidx,yidx]=ind2sub(size(zx),idx);
%     [xidx,yidx,zidx]=ind2sub(size(zx),idx);
    rad_tmp=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2);
%     rad_tmp=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2 + (zidx-refidx(3)).^2);
    radvec=[radvec;repmat(rad_tmp,acids_per_px(i),1)];
end
idx=find(zx>0);
[xidx,yidx]=ind2sub(size(zx),idx);
% data91_2.radius_acids_global=sqrt((xidx-refidx(1)).^2+(yidx-refidx(2)).^2);
data91_2.radius_acids_global=radvec;

figure;
imagesc(sum(data30.acidimg_global./1024,3));colorbar;
title('30 eV; ZX');
set(gca,'fontsize',30,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
caxis([0 0.12]);

figure;
imagesc(sum(data50.acidimg_global./1024,3));colorbar;
title('50 eV; ZX');
set(gca,'fontsize',30);
set(gcf,'Position',gcf_pos);
caxis([0 0.12]);

figure;
imagesc(reshape(sum(data30.acidimg_global./1024,1),50,50));colorbar;
title('30 eV; YX');
set(gca,'fontsize',30);
set(gcf,'Position',gcf_pos);
caxis([0 0.12]);

figure;
imagesc(reshape(sum(data50.acidimg_global./1024,1),50,50));colorbar;
title('50 eV; YX');
set(gca,'fontsize',30);
set(gcf,'Position',gcf_pos);
caxis([0 0.12]);


% figure;
% imagesc(sum(data80.acidimg_global,3));colorbar;

% figure;
% imagesc(sum(data91.acidimg_global,3));colorbar;

% figure;
% stem(x80,n80,'bd','linewidth',3.0,'MarkerSize',12);
% hold on;
% stem(x91,n91,'kd','linewidth',3.0,'MarkerSize',12);
% stem(x50,n50,'rd','linewidth',3.0,'MarkerSize',12);

% figure;

%% distribtions of acid x,y,z
clc;
close all;

idxtmp=find(data91.acidimg_global>0);
[zidx,xidx,yidx]=ind2sub(size(data91.acidimg_global),idxtmp);

figure;
hist(zidx);title('zidx')

figure;
hist(xidx);title('xidx');

figure;
hist(yidx);title('yidx');

%% Acid-generation radii

clc;
close all;
% 
dx=0.5;
vec=data30.radius_acids_global;
nbins=(max(vec)-min(vec))/dx;
[n30,x30]=hist(vec,nbins);
% x30=[0 x30];n30=[0 n30];
% n30=n30./sum(n30).*meanacids(1);
n30=n30./trapz(x30,n30);
figure;
hist(data30.radius_acids_global,25);
bar(x30,n30,1);
title(sprintf('30 eV; Mean = %.4f; Stdev = %.4f',mean(data30.radius_acids_global),std(data30.radius_acids_global)));
set(gca,'fontsize',30,'linewidth',3.0);
xlim([0 10]);
% ylim([0 0.4]);

vec=data50.radius_acids_global;
nbins=(max(vec)-min(vec))/dx;
[n50,x50]=hist(vec,nbins,nbins);
% x50=[0 x50];n50=[0 n50];
% n50=n50./sum(n50).*meanacids(2);
n50=n50./trapz(x50,n50);
figure;hist(data50.radius_acids_global,25);
bar(x50,n50,1);
title(sprintf('50 eV; Mean = %.4f; Stdev = %.4f',mean(data50.radius_acids_global),std(data50.radius_acids_global)));
set(gca,'fontsize',30,'linewidth',3.0);
xlim([0 10]);
% ylim([0 0.4]);

vec=data80.radius_acids_global;
nbins=(max(vec)-min(vec))/dx;
[n80,x80]=hist(vec,nbins);
% x80=[0 x80];n80=[0 n80];
% n80=n80./sum(n80).*meanacids(3);
n80=n80./trapz(x80,n80);
figure;hist(data80.radius_acids_global,25);
bar(x80,n80,1);
title(sprintf('80 eV; Mean = %.4f; Stdev = %.4f',mean(data80.radius_acids_global),std(data80.radius_acids_global)));
set(gca,'fontsize',30,'linewidth',3.0);
xlim([0 10]);
% ylim([0 0.4]);

vec=data91.radius_acids_global;
nbins=(max(vec)-min(vec))/dx;
[n91,x91]=hist(vec,nbins);
% x91=[0 x91];n91=[0 n91];
% n91=n91./sum(n91).*meanacids(4);
n91=n91./trapz(x91,n91);
figure;hist(data91.radius_acids_global,25);
bar(x91,n91,1);
title(sprintf('91 eV; Mean = %.4f; Stdev = %.4f',mean(data91.radius_acids_global),std(data91.radius_acids_global)));
set(gca,'fontsize',30,'linewidth',3.0);
xlim([0 10]);
% ylim([0 0.4]);

gcf_pos=[ 268   345   750   559];
gca_pos=[ 0.1300    0.1100    0.7750    0.8150];

figure;
plot(x30,n30./30,'-ob','linewidth',3.0);
hold on;
plot(x50,n50./50,'-ok','linewidth',3.0);
plot(x80,n80./80,'-or','linewidth',3.0);
plot(x91,n91./91,'-og','linewidth',3.0);
set(gca,'Position',gca_pos);
set(gca,'fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);
xlim([0 15]);

figure;
plot(x30,cumsum(n30),'-ob','linewidth',3.0);
hold on;
plot(x50,cumsum(n50),'-ok','linewidth',3.0);
plot(x80,cumsum(n80),'-or','linewidth',3.0);
plot(x91,cumsum(n91),'-og','linewidth',3.0);
set(gca,'Position',gca_pos);
set(gca,'fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);
xlim([0 15]);

%% Acid-generation radii [the polym, acid combined scenario]

clc;
% close all;
% 
dx=0.5;
vec=data30.radius_acids_global;
nbins=(max(vec)-min(vec))/dx;
[n30,x30]=hist(vec,nbins);
% x30=[0 x30];n30=[0 n30];
% n30=n30./sum(n30).*data30_2.meanacids(1);
n30=n30./sum(n30).*meanacids(1);
% n30=n30./sum(n30);
figure;
hist(data30.radius_acids_global,25);
bar(x30,n30./sum(n30),1);
title(sprintf('30 eV; Mean = %.4f; Stdev = %.4f',mean(data30.radius_acids_global),std(data30.radius_acids_global)));
set(gca,'fontsize',30,'linewidth',3.0);
xlim([0 10]);
ylim([0 0.32]);

vec=data50.radius_acids_global;
nbins=(max(vec)-min(vec))/dx;
[n50,x50]=hist(vec,nbins,nbins);
% x50=[0 x50];n50=[0 n50];
% n50=n50./sum(n50).*data50_2.meanacids(1);
n50=n50./sum(n50).*meanacids(2);
% n50=n50./sum(n50);
figure;hist(data50.radius_acids_global,25);
bar(x50,n50./sum(n50),1);
title(sprintf('50 eV; Mean = %.4f; Stdev = %.4f',mean(data50.radius_acids_global),std(data50.radius_acids_global)));
set(gca,'fontsize',30,'linewidth',3.0);
xlim([0 10]);
ylim([0 0.32]);

vec=data80.radius_acids_global;
nbins=(max(vec)-min(vec))/dx;
[n80,x80]=hist(vec,nbins);
% x80=[0 x80];n80=[0 n80];
% n80=n80./sum(n80).*data80_2.meanacids(1);
n80=n80./sum(n80).*meanacids(3);
% n80=n80./sum(n80);
figure;hist(data80.radius_acids_global,25);
bar(x80,n80./sum(n80),1);
title(sprintf('80 eV; Mean = %.4f; Stdev = %.4f',mean(data80.radius_acids_global),std(data80.radius_acids_global)));
set(gca,'fontsize',30,'linewidth',3.0);
xlim([0 10]);
ylim([0 0.32]);

vec=data91.radius_acids_global;
nbins=(max(vec)-min(vec))/dx;
[n91,x91]=hist(vec,nbins);
% x91=[0 x91];n91=[0 n91];
% n91=n91./sum(n91).*data91_2.meanacids(1);
n91=n91./sum(n91).*meanacids(4);
% n91=n91./sum(n91);
figure;hist(data91.radius_acids_global,25);
bar(x91,n91./sum(n91),1);
title(sprintf('91 eV; Mean = %.4f; Stdev = %.4f',mean(data91.radius_acids_global),std(data91.radius_acids_global)));
set(gca,'fontsize',30,'linewidth',3.0);
xlim([0 10]);
ylim([0 0.32]);

gcf_pos=[ 268   345   750   559];
gca_pos=[ 0.1300    0.1100    0.7750    0.8150];

figure;
plot(x30,n30,'-ob','linewidth',3.0);
hold on;
plot(x50,n50,'-ok','linewidth',3.0);
plot(x80,n80,'-or','linewidth',3.0);
plot(x91,n91,'-og','linewidth',3.0);
set(gca,'Position',gca_pos);
set(gca,'fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);
xlim([0 15]);
ylim([0 0.12]);

figure;
plot(x30,cumsum(n30)./1,'-ob','linewidth',3.0);
hold on;
plot(x50,cumsum(n50)./1,'-ok','linewidth',3.0);
plot(x80,cumsum(n80)./1,'-or','linewidth',3.0);
plot(x91,cumsum(n91)./1,'-og','linewidth',3.0);
set(gca,'Position',gca_pos);
set(gca,'fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);
xlim([0 15]);

%% Fitting the radius probability density function
% clc;
close all;

xvec=[x30 x50 x80 x91];
yvec=[n30 n50 n80 n91];

xvec=x91;
yvec=n91; 

[xvec,idx]=sort(xvec);
yvec=yvec(idx);
% yvec=yvec./sum(yvec);
yvec=yvec./trapz(xvec,yvec);

%%% fit an analytical function
meanval=2;
sigma=1;
radius_fit=0:0.1:20;
psf1=1/(sigma*sqrt(2*pi)).*exp(-(radius_fit-meanval).^2./(2*sigma^2));                  %psf1=psf1./max(psf1).*0.18;
psf2=1/(sigma*sqrt(2*pi)).*radius_fit.*exp(-(radius_fit-meanval).^2./(2*sigma^2));      %psf2=psf2./max(psf2).*0.18;

%%% max,meanval,sigma:
Rayl=@(parms,radius) 1/(parms(3)*sqrt(2*pi)).*radius.*parms(1).*exp(-(radius-parms(2)).^2./(2*parms(3)^2));
% Rayl2=@(parms,radius) 1/(parms(2)*sqrt(2*pi)).*radius.*exp(-(radius-parms(1)).^2./(2*parms(2)^2));
% Rayl2=@(parms,radius) 1/(parms(2))^2.*radius.*exp(-(radius-parms(1)).^2./(2*parms(2)^2));
% Rayl2=@(parms,radius) 1/(parms(1))^2.*radius.*exp(-(radius).^2./(2*parms(1)^2));
% Rayl2=@(parms,radius) parms(2)/(parms(1))^2.*radius.*exp(-(radius).^2./(2*parms(1)^2));
Rayl2=@(parms,radius) 1/(parms(1))^2.*radius.*exp(-(radius).^2./(2*parms(1)^2));

Gauss=@(parms,radius) 1/(parms(2)*sqrt(2*pi)).*exp(-(radius-parms(1)).^2./(2*parms(2)^2));

maxval=0.03;
meanval=1;
sigma=1;

parms=[meanval sigma  ];
% parms=[meanval sigma];
% parms=[sigma maxval];
parms=[sigma];

options=optimset('lsqcurvefit');
options.MaxFunEvals=1000;
% parms=lsqcurvefit(@(parms,radius) 1/(parms(3)*sqrt(2*pi)).*radius.*parms(1).*exp(-(radius-parms(2)).^2./(2*parms(3)^2)),parms,xvec,yvec,[],[],options);
% parms=lsqcurvefit(@(parms,radius) 1/(parms(2)*sqrt(2*pi)).*radius.*exp(-(radius-parms(1)).^2./(2*parms(2)^2)),parms,xvec,yvec,[],[],options);
% parms=lsqcurvefit(@(parms,radius) 1/(parms(2))^2.*radius.*exp(-(radius-parms(1)).^2./(2*parms(2)^2)),parms,xvec,yvec,[],[],options);
% parms=lsqcurvefit(@(parms,radius) parms(2)/(parms(1))^2.*radius.*exp(-(radius).^2./(2*parms(1)^2)),parms,xvec,yvec,[],[],options);
parms=lsqcurvefit(@(parms,radius) 1/(parms(1))^2.*radius.*exp(-(radius).^2./(2*parms(1)^2)),parms,xvec,yvec,[],[],options);
% parms=lsqcurvefit(@(parms,radius) 1/(parms(2)*sqrt(2*pi)).*exp(-(radius-parms(1)).^2./(2*parms(2)^2)),parms,xvec,yvec,[],[],options);

gcf_pos=[680   331   817   647];
gca_pos=[0.1300    0.1100    0.7750    0.8150];

mse=sqrt(mean((Rayl2(parms,xvec)-yvec).^2))
% mse=sqrt(mean((Gauss(parms,xvec)-yvec).^2))

figure;
plot(xvec,yvec,'ob','MarkerSize',12,'linewidth',2.0,'linewidth',3.0)
hold on;
plot(radius_fit,Rayl2(parms,radius_fit),'--k','linewidth',4.0);
% plot(radius_fit,Gauss(parms,radius_fit),'--k','linewidth',4.0);
title(sprintf('Sigma = %.4f',parms(1)));
% title(sprintf('Mean = %.4f; Sigma = %.4f',parms(1),parms(2)));
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
% ylim([0 0.5]);
xlim([0 10]);
