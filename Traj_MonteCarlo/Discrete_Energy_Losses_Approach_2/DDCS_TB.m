%% DDCS Testbench
clear;
% clc;
% close all;

addpath('F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\GlobalFunctions');

na=1.2/120*6.02*1e23;
JeV=1.6*1e-19;
E=[10 20 30 40 50 60 70 80 90 100 200];
E=[5 10 15 20 30 40 50 60 70 80 90 100 200 1000];
E=[6 11 15 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000];
E=[6:100 150:50:1000];
E=5:1:200;
E=[6:200 225:25:1000];
E=[16:5:710];
E=[50];

ddcs_sum2=[];
ics=[];

fig1=figure;

Elossmat=[];
dsigdE=[];
dsigdOmega=[];

tstart_global=tic;
for i = 1:length(E)
    tstart=tic;
    fprintf('Energy %d of %d\n',i,length(E));
%     [ddcs,dsigdOmega(i,:),ics(i,:),sigTot(i)]=ddcscalc(E(i));
    outparms=ddcscalc(E(i));
%     outparms2=dsigdEcalc(E(i));
%     Elossmat_2(i,:)=outparms2.Elossmat;
%     dsigdE_2(i,:)=outparms2.dsigdE;
%     imfp_2(i,:)=outparms2.imfp;
    
    ELF=outparms.ELF;
    qvec=outparms.qvec;
    dsigdOmega_1(i,:)=outparms.dsigdOmega;
    ics(i,:)=outparms.ics;
    sigTot(i)=outparms.sigTot;
    Sp(i)=outparms.Sp;
%     Sp2(i)=outparms.Sp2;
    imfp(i)=outparms.imfp;
    Eloss=outparms.Eloss;
    theta=outparms.theta;
    theta_d=theta.*180/pi;
        
    ddcs=outparms.ddcs; %unit: m^2/J
%     ddcs=ddcs.*1e4*JeV;  % converted to cm2/eV;
    theta2=reshape(theta,[size(ddcs,1),1]);
    theta2=repmat(theta2,[1 size(ddcs,2)]);
    d2sigdTdw=2*pi.*sin(theta2).*ddcs;
        
    for theta_count=1:length(theta)
        thetamat(i,theta_count)=theta(theta_count);
        dsigdOmega(i,theta_count)=trapz(Eloss,ddcs(theta_count,:)); % DDCS is /m/eV/steradian
    end
    
    for Eloss_count=1:length(Eloss)
        Elossmat(i,Eloss_count)=Eloss(Eloss_count);
        dsigdE(i,Eloss_count)=trapz(theta',2*pi.*sin(theta').*ddcs(:,Eloss_count)); % DDCS is /m/eV/steradian
    end
%     dsigdE(i,:)=dsigdE(i,:)./max(dsigdE(i,:));
    
%     sigTot2(i)=trapz(Eloss,dsigdE(1,:));
    
%     icsdata{i}.E=outparms.CrossSect.E;
%     icsdata{i}.theta=outparms.CrossSect.theta;
%     icsdata{i}.CS_2D=outparms.CrossSect.CS_2D;
    
    ddcs_sum2(i)=sum(sum(ddcs));
    
    figure(fig1);
%     figure;
%     imagesc(ddcs);
%     imagesc(Eloss,theta_d,d2sigdTdw);
    imagesc(Eloss,theta_d,d2sigdTdw);
    ylabel('theta (degrees)');
    xlabel('Eloss (eV)');
    title(['E = ' num2str(E(i)) ' eV']);
    colorbar;
    drawnow;
    
    tloop(i)=toc(tstart);
    fprintf('... Took %.2f s\n',tloop(i));
end
tend_global=toc(tstart_global);
fprintf('\nTotal Time = %.4f s = %.4f min\n',tend_global,tend_global/60);

%%%% The Sp data from Truong (2015);
truong_data_path='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\Traj_MonteCarlo\Discrete_Energy_Losses_Approach\';
truong_data=load([truong_data_path 'TruongData_PMMA_PS.mat']);

figure;plot(E,Sp,'linewidth',3.0);
if exist('truong_data')
    hold on;
    plot(truong_data.PMMA.E,truong_data.PMMA.Sp,'--ob');
    plot(truong_data.PS.E,truong_data.PS.Sp,'--or');
    legend('PHS-mermin','Truong-PMMA','Truong-PS');
end

xlabel('Energy (eV)','fontsize',20);
ylabel('Normalized Sp','fontsize',20);
set(gca,'XScale','linear');
xlim([0 1000]);
set(gca,'fontsize',20,'linewidth',3.0);

%% relative probabilities over various energy regimes
% clc;
close all;

Etmp=Elossmat(1,:);
pdf=dsigdE(1,:);

int=[];
int(1)=trapz(Etmp(Etmp>0 & Etmp<=5),pdf(Etmp>=0 & Etmp<=5));
int(2)=trapz(Etmp(Etmp>5 & Etmp<=12),pdf(Etmp>=5 & Etmp<=12));
int(3)=trapz(Etmp(Etmp>12 & Etmp<=max(Etmp)),pdf(Etmp>=12 & Etmp<=max(Etmp)));

int=[int.*1e-9;int./sum(int)];

%% Run this to save the data
Ef=15.5;

% save('DDCSData\DDCSdata_Fuji_Ef=15.5_Elossmin=3eV_Erange=[19,200].mat','Elossmat','dsigdE','thetamat','dsigdOmega','E','Sp');
save('DDCSData\DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,710].mat','Elossmat','dsigdE','thetamat','dsigdOmega','E','Sp');
% save('DDCSData\DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,200].mat','Elossmat','dsigdE','thetamat','dsigdOmega','E','Sp','Elossmat_2','dsigdE_2');

% save('DDCSData\Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=3eV_Erange=[19,200]_DDCSData.mat','E','imfp','Sp','Ef');
save('DDCSData\Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,710]_DDCSData.mat','E','imfp','Sp','Ef');
% save('DDCSData\Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData.mat','E','imfp','Sp','Ef','imfp_2');

%% Test the ddcs through video
clc;
close all;

fig1=figure;
ics=[];
for i = 1:size(dsigdOmega,1)
    ics(i,:)=cumsum(2*pi.*sin(theta).*dsigdOmega(i,:)).*(theta(2)-theta(1));
end

