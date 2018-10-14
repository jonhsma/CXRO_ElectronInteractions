%% Scattering Sim: Get the stopping power
clear;
% clc;
close all;

addpath('F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\GlobalFunctions');

cross_sect_data_path='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\Traj_MonteCarlo\CrossSect_Data\';

scattdata.vibr=load([cross_sect_data_path 'VibrExcit_Data_Khakoo_2013.mat']);
scattdata.vibr.datasrc='Khakoo';

% scattdata.vibr.datasrc='Frohlich';
% % scattdata.vibr.E=0:0.1:500;
% scattdata.vibr.eps0=1.5^2;
% scattdata.vibr.epsinf=1;
% scattdata.vibr.hbarw=[0.1 0.1];

% scattdata.vibr.imfp=eph_imfp(scattdata.vibr);
scattdata.vibr.imfp_func=@ephscatt;

% optdata_path='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\IMFP_Analysis\MatFiles\IMFP_Components\';
optdata_path='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';
% scattdata.optical=load([optdata_path 'Sp_IMFP_Inelastic_Components_Ef=0p5eV_Elossmin=0.001eV_Erange=[5,1000]_DDCSData.mat']); % _v2 has data structures better suited to the subsequently called programs
% scattdata.optical=load([optdata_path 'Sp_withICSData_IMFP_Inelastic_Components_Ef=10eV_Elossmin=0.001eV_Erange=[5,1000]_DDCSData.mat']);
% scattdata.optical=load([optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData_EQCons=Pines.mat']);
scattdata.optical=load([optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData.mat']);

%%% File containing the Sp values
% Spdata_path='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\Traj_MonteCarlo\Discrete_Energy_Losses_Approach\';
% Spdata_file='Sp_IMFP_Inelastic_Components_Ef=5eV_Elossmin=0.001eV_Erange=[6,200]_DDCSData.mat';
% Spdata=load([Spdata_path Spdata_file]);

%%% File containing the imfp values
% imfpdata_path='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\IMFP_Analysis\';
% imfpdata_file='IMFPData_Ef=10eV_Erange=[5,1000]_Elossmin=0.001eV.mat';
% imfpdata=load([imfpdata_path imfpdata_file]);

% scattdata.optical.E=Spdata.E;
% scattdata.optical.Sp=Spdata.Sp;
% scattdata.optical.imfp=imfpdata.imfp;
% scattdata.optical.imfp=interp1(imfpdata.E,imfpdata.imfp,scattdata.optical.E);
% scattdata.optical.Ef=imfpdata.Ef;

TruongData=load('TruongData_PMMA_PS.mat');

% inel_dcsdata_path='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\Traj_MonteCarlo\';
% scattdata.optical.inel_dcsdata=load([inel_dcsdata_path 'Inelastic_DCS_Data.mat']);
% [filename,pathname,fid]=uigetfile(['*.mat'],'Select DDCS Data File');
pathname='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';
filename='DDCSdata_Ef=0p5_Elossmin=0.001eV_Erange=[5,1000].mat';
filename='DDCSdata_withICSData_Ef=10_Elossmin=0.001eV_Erange=[5,1000].mat';
filename='DDCSdata_Fuji_Ef=15.5_Elossmin=3eV_Erange=[19,200]_EQCons=Pines.mat';
filename='DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,200].mat';
scattdata.optical.inel_dcsdata=load([pathname filename]);

scattdata.E_inel_thr=min(scattdata.optical.E);

%%%%%% create PAG grid
rho_pag=0.4; % per nm3
rho_polym=6; % / nm3
% rho_pag=0;
px_nm=[1 1 1]*1;

univ.size_nm=[50 50 50]; % nm
univ.px_nm=px_nm;
univ.npx=round(univ.size_nm./univ.px_nm);

xgrid=-univ.size_nm(1)/2:px_nm(1):univ.size_nm(1)/2-px_nm(1);
ygrid=-univ.size_nm(2)/2:px_nm(2):univ.size_nm(2)/2-px_nm(2);
zgrid=-univ.size_nm(3)/2:px_nm(3):univ.size_nm(3)/2-px_nm(3);
[xgrid,zgrid,ygrid]=meshgrid(xgrid,xgrid,zgrid);
univ.grid.x=xgrid;
univ.grid.y=ygrid;
univ.grid.z=zgrid;

pagimg=rho_pag*prod(univ.px_nm).*ones(univ.npx);
pagimg=poissrnd(pagimg);
elecimg=zeros(size(pagimg));

%%%%%% Initialize the event structure [holds info about excitations triggered]
% event{1}.xyz=round(univ.npx/2);
event{1}.xyz=[0 0 0];
xyzglobal.x=event{1}.xyz(1);
xyzglobal.y=event{1}.xyz(2);
xyzglobal.z=event{1}.xyz(3);

E_in=91;

event{1}.Ein=E_in;
event{1}.Eout=0;

event{1}.Eloss=0;
event{1}.Ese=91;
event{1}.act='SE';
event{1}.nacid=0;
event{1}.nSE=0;
event{1}.elecimg=elecimg;
event{1}.pag.img=pagimg;
event{1}.pag.rcnrad=2; % nm
event{1}.pag.rho=rho_pag;
% event{1}.controlparms.model_pag_devel=1;
event{1}.univ=univ;
% event{1}.scatt_Elim=20;
% event{1}.lowEthr=0;
% event{1}.lowEimfp=3.5;

% ntrials=sum(absimg(:));
% abspos=find(absimg~=0);
ntrials=1024;

plotfig=0;
if plotfig~=0
    fig1=figure;hold on;colors='brgk';
end

tstart=tic;

Esweep=linspace(20,100,5);
Esweep=[29];

pathlen=[];
Energy=[];

%%%% electron dose parameters
elecimg_inc=zeros([univ.npx(1) univ.npx(2)]);
pattsize=5; % nm
rowsel=round(univ.npx(1)/2-pattsize/2/univ.px_nm(1)):round(univ.npx(1)/2+pattsize/2/univ.px_nm(1));
colsel=round(univ.npx(2)/2-pattsize/2/univ.px_nm(2)):round(univ.npx(2)/2+pattsize/2/univ.px_nm(2));
% rowsel=rowsel(2:end); % remove the extra 1 addeded due to integer math
% colsel=colsel(2:end);
elecimg_inc(rowsel,colsel)=1;
elecimg_bin=elecimg_inc;

% [78,157,314];
Dose=4; % /pixel;

elecimg_inc=Dose.*elecimg_inc;
% elecimg_inc=floor(elecimg_inc);
% rng(56);
rng('shuffle');
% elecimg_inc=poissrnd(elecimg_inc);
rng('shuffle');

%%% if you want to do a trial sweep and have only 1 electron incident at
%%% center, uncomment the 2 lines below:
%     elecimg_inc=zeros([univ.npx(1) univ.npx(2)]);
%     elecimg_inc(round(size(elecimg_inc,1)/2),round(size(elecimg_inc,2)/2))=Dose;

figure;imagesc(elecimg_inc);colorbar;title('Electron image');drawnow;

[x_inc,y_inc]=ind2sub(size(elecimg_inc),find(elecimg_inc>=1));

logfile_fid=fopen('logfile.dat','w');
fprintf(logfile_fid,'Simulation started at %s\n',datestr(now));
for E_count=1:length(Esweep)
    event{1}.Ein=Esweep(E_count);
    event{1}.Ese=Esweep(E_count);
        
    eventdata_global={};
    Eloss=[];

    nacids_total=[];
    nacids_total2=[];
    Energy_global=[];
    imfp_global=[];
    steplen_global=[];
    radius_acids_global=[];
    radius_ions_global=[];
    act_global={};
    theta_global=[];
    phi_global=[];
    xyz_electron_global=[];
    acidimg_global=zeros(univ.npx);
    ionimg_global=zeros(univ.npx);
    radius_acids_thrutrial={};
    radius_ions_thrutrial={};
    nacids_thrutrial=[];
    nions_thrutrial=[];
    meanacids_thrutrial=[];
    meanions_thrutrial=[];
    
    theta_init=[];

    for trial_count=1:ntrials
        acid_xyz=[];
        radius_acids=[];
        radius_ions=[];
        ion_xyz=[];

        tst_trials=tic;
        fprintf(logfile_fid,'\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
        fprintf('\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
        
        pagimg=rho_pag*prod(univ.px_nm).*ones(univ.npx);
        polym_img=rho_polym*prod(univ.px_nm).*ones(univ.npx);
        
        rng('shuffle');
%         rng(56);
        pagimg=poissrnd(pagimg);
        
        rng('shuffle');
        polym_img=poissrnd(polym_img);
        rng('shuffle');
        
        pagimg_orig=pagimg;
        polymimg_orig=polym_img;
        fprintf(logfile_fid,'Initial <PAG> = %.4f/nm^3\n',mean(pagimg_orig(:))/prod(univ.px_nm));
        fprintf(logfile_fid,'Initial <Polymer> = %.4f/nm^3\n',mean(polymimg_orig(:))/prod(univ.px_nm));
        fprintf('\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
%         event{1}.pag.img=pagimg;
        pagdata.pagimg=pagimg;
        pagdata.acid_act_xyz_idx=[];
        polymdata.polym_img=polym_img;
        polymdata.SE_act_xyz_idx=[];
        
        acidimg=zeros(size(pagimg));
        
%         xyz_acids=[];radius_acids=[];
        xyz_electron=[];
%         nacids_thrutrial=[];

        for xyz_count=1:length(x_inc)
            tst_coordinates=tic;
            fprintf(logfile_fid,'\n...Co-ordinate %d of %d\n',xyz_count,length(x_inc));
%             fprintf('\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
            nelectrons=elecimg_inc(x_inc(xyz_count),y_inc(xyz_count));
            xval=univ.grid.x(1,x_inc(xyz_count),y_inc(xyz_count));
            yval=univ.grid.y(1,x_inc(xyz_count),y_inc(xyz_count));
            zval=min(univ.grid.z(:)); % put the electron at the origin in z.
%             zval=0; % put it in the center

            for photon_count=1:nelectrons
                tst_nelec=tic;
                fprintf('......Electron %d of %d\n',photon_count,nelectrons);
                fprintf(logfile_fid,'......Electron %d of %d\n',photon_count,nelectrons);
%                 fprintf('\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
                event{1}.xyz=[xval yval zval]; % [z,x,y]
                
                %%%%% initial conditions [for LEEM, set both =0]:
%                 event{1}.theta_in=0;
%                 event{1}.phi_in=0;
                
                %%%%% initial conditions [for non-LEEM: do whatever you want, perhaps photo-electron emission angle]:
                rng('shuffle');
%                 event{1}.theta_in=-pi+pi*rand;
                event{1}.theta_in=0;
                theta_init=[theta_init event{1}.theta_in];
                rng('shuffle');
                event{1}.phi_in=2*pi*rand;
                
%                 event{1}.pag.img=pagimg;
                
                xyzglobal.x=event{1}.xyz(1);
                xyzglobal.y=event{1}.xyz(2);
                xyzglobal.z=event{1}.xyz(3);
%                 xyzglobal.z=0; % put the electron at the center.

                eventdata={};
                pagimg_pre=pagdata.pagimg;
                polym_img=polymdata.polym_img;
                [eventdata,pagdata,polymdata]=Scattcalc_lowE(event,scattdata,eventdata,xyzglobal,pagdata,polymdata,logfile_fid);
    %             eventdata_global{trial_count,E_count}=eventdata;

%                 acidimg=acidimg+pagimg-eventdata{end}.pag.img;
                acidimg=acidimg+pagimg_pre-pagdata.pagimg;
                if any(pagimg_pre(:)<pagdata.pagimg(:))
                    dbg=1;
                end
%                 pagimg=eventdata{end}.pag.img;
%                 pagimg(pagimg<0)=0;
                nacids_total=[nacids_total sum(acidimg(:))];

                act={};nSE=[];imfp=[];theta=[];phi=[];nacids=[];

                pathlen(trial_count,1)=0; % its 0 when electron energy is the incident energy.
                nacids_count=1;
%                 xyz_acids=[];%radius_acids=[];
                for count = 1:length(eventdata)
                    Ein(trial_count,count)=eventdata{count}.Ein;
                    Eloss(trial_count,count)=eventdata{count}.Eloss;
                    Ese(trial_count,count)=eventdata{count}.Ese;
%                     xyz(count,:)=eventdata{count}.xyz;
                    xyz_electron=[xyz_electron;eventdata{count}.xyz];
                    act{count}=eventdata{count}.act;
                    nSE(count)=eventdata{count}.nSE;
                    nacids(count)=eventdata{count}.nacid;
                    pathlen(trial_count,count+1)=eventdata{count}.pathlen;
                    Energy(trial_count,count)=eventdata{count}.Ein;
        %             Energy=[Energy eventdata{end}.Eout];
        
                    imfp(count)=eventdata{count}.imfp;
                    theta(count)=eventdata{count}.theta;
                    phi(count)=eventdata{count}.phi;
                    
                    Energy_global=[Energy_global;eventdata{count}.Ein];
                    imfp_global=[imfp_global;eventdata{count}.imfp];
                    steplen_global=[steplen_global;eventdata{count}.rnew];
                    act_global{length(act_global)+1}=eventdata{count}.act;
                    theta_global(length(theta_global)+1)=eventdata{count}.theta;
                    phi_global(length(phi_global)+1)=eventdata{count}.phi;
                    xyz_electron_global(size(xyz_electron_global,1)+1,1:3)=eventdata{count}.xyz;
                    
                    if plotfig~=0
                        figure(fig1);
                        switch act{count}
                            case 'SE'
                                plot3(eventdata{count}.xyz(1),eventdata{count}.xyz(2),eventdata{count}.xyz(3),'ob','markersize',12,'linewidth',3.0);
                            case '6eVRes'
                                plot3(eventdata{count}.xyz(1),eventdata{count}.xyz(2),eventdata{count}.xyz(3),'or','markersize',12,'linewidth',3.0);
                            case 'vibr'
                                plot3(eventdata{count}.xyz(1),eventdata{count}.xyz(2),eventdata{count}.xyz(3),'ok','markersize',12,'linewidth',3.0);
                            case 'none'
                                plot3(eventdata{count}.xyz(1),eventdata{count}.xyz(2),eventdata{count}.xyz(3),'og','markersize',12,'linewidth',3.0);
                            case '6eVRes-none'
                                plot3(eventdata{count}.xyz(1),eventdata{count}.xyz(2),eventdata{count}.xyz(3),'*r','markersize',12,'linewidth',3.0);
                        end
                        drawnow;
%                         hold on;plot3(xyz(:,1),xyz(:,2),xyz(:,3),'linewidth',3.0)
                    end
                end
                Energy(trial_count,count+1)=eventdata{count}.Eout;
                nacids_total2=[nacids_total2 sum(nacids)];
                
                tend=toc(tst_nelec);
                fprintf(logfile_fid,'......Took %.2f s \n',tend);
            end
            
            tend_coordinates=toc(tst_coordinates);
            fprintf(logfile_fid,'...Took %.2f s\n',tend_coordinates);
        end
        
        %%% Do acid counting etc. right before next trial and energy values
        %%% are used:
        acid_act_xyz_idx=pagdata.acid_act_xyz_idx;
        acidimg_global(acid_act_xyz_idx)=acidimg_global(acid_act_xyz_idx)+1;
        for acid_count=1:length(acid_act_xyz_idx)
            [actidx_x,actidx_y,actidx_z]=ind2sub(size(pagdata.pagimg),acid_act_xyz_idx(acid_count));
            xyz_tmp=[univ.grid.x(actidx_x,actidx_y,actidx_z) univ.grid.y(actidx_x,actidx_y,actidx_z) univ.grid.z(actidx_x,actidx_y,actidx_z)];
            acid_xyz=[acid_xyz;xyz_tmp];
        end
        
        SE_act_xyz_idx=polymdata.SE_act_xyz_idx;
        ionimg_global(SE_act_xyz_idx)=ionimg_global(SE_act_xyz_idx)+1;
        for ion_count=1:length(SE_act_xyz_idx)
            [ionidx_x,ionidx_y,ionidx_z]=ind2sub(size(polymdata.polym_img),SE_act_xyz_idx(ion_count));
            xyz_tmp=[univ.grid.x(ionidx_x,ionidx_y,ionidx_z) univ.grid.y(ionidx_x,ionidx_y,ionidx_z) univ.grid.z(ionidx_x,ionidx_y,ionidx_z)];
            ion_xyz=[ion_xyz;xyz_tmp];
        end
        
        meanions_thrutrial(trial_count)=mean(length(SE_act_xyz_idx));
        nions_thrutrial(trial_count)=length(SE_act_xyz_idx);
        if nions_thrutrial(trial_count)>0
            radius_ions=[radius_ions;sqrt((ion_xyz(:,1)-xval).^2+(ion_xyz(:,2)-yval).^2+(ion_xyz(:,3)-zval).^2)];
            radius_ions_thrutrial{trial_count}=radius_ions;
            radius_ions_global=[radius_ions_global;radius_ions];
        else
            radius_ions_thrutrial{trial_count}=NaN;
        end
        
        meanacids_thrutrial(trial_count)=mean(length(acid_act_xyz_idx));
        nacids_thrutrial(trial_count)=length(acid_act_xyz_idx);
        if nacids_thrutrial(trial_count)>0
            radius_acids=[radius_acids;sqrt((acid_xyz(:,1)-xval).^2+(acid_xyz(:,2)-yval).^2+(acid_xyz(:,3)-zval).^2)];
            radius_acids_thrutrial{trial_count}=radius_acids;
            radius_acids_global=[radius_acids_global;radius_acids];
        else
            radius_acids_thrutrial{trial_count}=NaN;
        end
        dbg=1;
        
%         save(sprintf('LEEM_SimResults_4_04302017\\Ein=%.2f_Dose=%.2f_Ef=15.5_pag-Emin=5_rcnrad=2_PAG=0.4_T%d.mat',Esweep(E_count),Dose,trial_count));
%         save(sprintf('LEEM_SimResults_low10AllPAGs\\Ein=%.2f_Dose=%.2f_Ef=15.5_pag-Emin=5_rcnrad=2_PAG=0.4_T%d.mat',Esweep(E_count),Dose,trial_count));
%         save(sprintf('Scratch_Folder\\Ein=%.2f_Dose=%.2f_Ef=15.5_pagEa=5_Scatt-Elim=%.2f_rcnrad=2_PAG=0.4_T%d.mat',Esweep(E_count),Dose,event{1}.scatt_Elim,trial_count));
        save(sprintf('Scratch_Folder\\Ein=%.2f_Dose=%.2f_Ef=15.5_pagEa=5_Scatt-Elim=%.2f_rcnrad=2_PAG=0.4_T%dB.mat',Esweep(E_count),Dose,20,trial_count));

%         tend=toc(tst_trials);
%         fprintf('...Took %.2f s \n',tend);
    end
    
%     save(sprintf('QE_Stats_Sims_7\\QEStats_Dose=%.2f_RcnRad=%.2f_PAG=%.2f_E=%.2feV_pagEa=5_Scatt-Elim=20_PAG-EaMin=5_%dTrials_RcnThenMove.mat',Dose,event{1}.pag.rcnrad,rho_pag,Esweep(E_count),ntrials));
%     save(sprintf('QE_Stats_low10AllPAGs\\QEStats_Dose=%.2f_RcnRad=%.2f_PAG=%.2f_E=%.2feV_pagEa=5_Scatt-Elim=20_PAG-EaMin=5_%dTrials_RcnThenMove.mat',Dose,event{1}.pag.rcnrad,rho_pag,Esweep(E_count),ntrials));
    %     save(sprintf('QEStats_Dose=%.2f_RcnRad=%.2f_PAG=%.2f_E=%.2feV_pagEa=5_Scatt-Elim=5_PAG-EaMin=5_%dTrials_ReactThenMove_RandomAngles.mat',Dose,event{1}.pag.rcnrad,rho_pag,Esweep(E_count),ntrials));

    tend=toc(tstart);
    fprintf(logfile_fid,'Total time = %.4f s = %.4f min = %.4f hours\n',tend,tend/60,tend/3600);
    fprintf('Total time = %.4f s = %.4f min = %.4f hours\n',tend,tend/60,tend/3600);

    fprintf(logfile_fid,'Total # of acids = %d\n',sum(nacids));
    fprintf('Total # of acids = %d\n',sum(nacids));

    fprintf(logfile_fid,'%d; <Acids> = %.4f; sig_Acids = %.4f\n',ntrials,mean(nacids_thrutrial),std(nacids_thrutrial));
    fprintf(1,'%d trials; <Acids> = %.4f; sig_Acids = %.4f\n',ntrials,mean(nacids_thrutrial),std(nacids_thrutrial));

end

%% Save QE stats files
clc;
close all;

% save(sprintf('QEStats_RcnRad_PAgLoading_Sweep\\QEStats_RcnRad=%.2f_PAG=%.2f_E=80eV_%dTrials.mat',event{1}.pag.rcnrad,rho_pag,ntrials));
% save(sprintf('QEStats_RcnRad_PAgLoading_Sweep\\QEStats_RcnRad=%.2f_PAG=%.2f_E=80eV_pagEa=3_%dTrials.mat',event{1}.pag.rcnrad,rho_pag,ntrials));
% save(sprintf('ThruDose_DeltaFunc\\QEStats_Dose=%.2f_RcnRad=%.2f_PAG=%.2f_E=80eV_pagEa=3_%dTrials.mat',Dose,event{1}.pag.rcnrad,rho_pag,ntrials));
save(sprintf('QEStats_Sims_2\\QEStats_Dose=%.2f_RcnRad=%.2f_PAG=%.2f_E=%.2feV_pagEa=3_Scatt-Elim=20_PAG-EaMin=5_%dTrials.mat',Dose,event{1}.pag.rcnrad,rho_pag,Esweep,ntrials));

%% Look through act_global to find certain events
clc;
close all;

idx=[];idxcount=1;
for i = 1:length(act_global)
    if strcmp(act_global{i},'6eVRes-none')==1
        idx(idxcount)=i;
        idxcount=idxcount+1;
    end
end

%% 2D matrix containing distance of electron reach for a given (x,y)
clc;
close all;

elecimg=zeros([univ.npx(2) univ.npx(3)]);
xmat=univ.grid.x;
ymat=univ.grid.y;
zmat=univ.grid.z;
tstart=tic;
fprintf('Creating electron surface\n');
for i = 1:size(xyz_electron,1)
    fprintf('Step %d of %d\n',i,size(xyz_electron,1));
    xtmp=xyz_electron(i,1);
    ytmp=xyz_electron(i,2);
    ztmp=xyz_electron(i,3);
    dist=(xmat-xtmp).^2+(ymat-ytmp).^2+(zmat-ztmp).^2;
    [tmpidx]=find(dist==min(dist(:)));
    [zidx,xidx,yidx]=ind2sub(size(xmat),tmpidx);
    xidx=xidx(1);yidx=yidx(1);zidx=zidx(1);
    elecimg(yidx,xidx)=xyz_electron(i,3);
end
tend=toc(tstart);
fprintf('Completed creating the electron image front, took %.4f s\n',tend);


% fprintf('Creating electron surface\n');
% for i = 1:size(xyz_electron,1)
%     fprintf('Step %d of %d\n',i,size(xyz_electron,1));
%     xtmp=xyz_electron(i,1);
%     ytmp=xyz_electron(i,2);
%     ztmp=xyz_electron(i,3);
%     dist=(xmat-xtmp).^2+(ymat-ytmp).^2+(zmat-ztmp).^2;
%     [tmpidx]=find(dist==min(dist(:)));
%     [zidx,xidx,yidx]=ind2sub(size(xmat),tmpidx);
%     xidx=xidx(1);yidx=yidx(1);zidx=zidx(1);
%     elecimg(yidx,xidx)=xyz_electron(i,3);
% end
% tend=toc(tstart);
% fprintf('Completed creating the electron image front, took %.4f s\n',tend);

%% Sp analysis from the calculated result
clc;
close all;

numtrials=size(Energy,1);
% numtrials=10;
dEdz=[];E=[];dEdz_int=[];cumul_path=[];pathlen_int=[];
cumul_path_int=[];
% Efit=1:1000; % energy values for fitting the dE/dz
Efit_min=5;
Efit_max=min(Energy(:,2));
Efit_max=900;
Efit=Efit_min:Efit_max;
good_count=1;
bad_idx=[];
for count=1:numtrials % sweep through all the trials
% for count=74
    fprintf('Calculating Monte Carlo Sp Trial %d of %d\n',count,numtrials);
    Evec=Energy(count,:);
    Evec=Evec(1:end-1);
    Elossvec=Eloss(count,:);
    pathvec=pathlen(count,2:end);
    idx=find(Elossvec~=0);
    Elossvec=Elossvec(idx);
    pathvec=pathvec(idx);
    Evec=Evec(idx);
    
    if length(pathvec)<=1
        continue;
    end
    
    %%%% the derivative approach [incorrectly, actually for Sp! Do
    %%%% forward-looking one below.
%     pathvec=cumsum(pathvec);
%     dEdzvec=abs(deriv(pathvec,Evec));
    
    %%%% Sp through delta_E/delta_path
    dEdzvec=Elossvec./pathvec;
    
%     dEdzvec=1./abs(deriv(Evec,pathvec)); % Alternative approach
    dEdz(good_count,1:length(dEdzvec))=dEdzvec;
    E(good_count,1:length(Evec))=Evec;
    cumul_path(good_count,1:length(pathvec))=cumsum(pathvec);
    cumul_path_int(good_count,:)=interp1(Evec,cumsum(pathvec),Efit,'linear');
    
    if min(Efit)<min(Evec)
        bad_idx=[bad_idx count];
    end
    
    dEdz_int(good_count,:)=interp1(Evec,abs(dEdzvec),Efit,'linear','extrap');
    pathlen_int(good_count,:)=interp1(Evec,cumsum(pathvec),Efit,'linear','extrap');
    good_count=good_count+1;
end

%%%%% Calculate the average dE/dz by binning
Ebin_size=0.5; % eV
Ebin_min=min(E(E~=0));
Ebin_max=max(E(:));
% numbins=100;
Efit_vec=Ebin_min+Ebin_size:2*Ebin_size:Ebin_max;
dEdz_mean=[];
for i = 1:length(Efit_vec)
    dEdz_tmp=dEdz(E>=Efit_vec(i)-Ebin_size & E<=Efit_vec(i)+Ebin_size);
    dEdz_mean(i)=mean(dEdz_tmp(~isnan(dEdz_tmp)));
end

cumul_path_int_mean=mean(cumul_path_int,1);
dEdz_2=abs(deriv(cumul_path_int_mean,Efit));

figure;
plot(cumsum(pathlen,2)',Energy','ob','MarkerSize',12,'linewidth',3.0)
hold on; 
plot(cumsum(pathlen,2)',Energy','k','linewidth',3.0)
xlabel('Cumul. path Length (nm)','fontsize',30);
ylabel('Energy (eV)','fontsize',30);
set(gca,'fontsize',30,'linewidth',3.0);

figure;
% plot(Efit',dEdz_int','-','linewidth',1.0);
% plot(Evec,abs(dEdzvec),'ob','linewidth',3.0);
plot(E(:,2:end)',dEdz(:,2:end)','ob','linewidth',3.0,'MarkerSize',12);
hold on;
% plot(Efit_vec,dEdz_mean,'-k','linewidth',3.0);
% plot(Efit,mean(dEdz_int,1),'-k','linewidth',3.0)
% plot(Efit',mean(dEdz_int,1),'-b','linewidth',6.0);
% errorbar(Efit',mean(dEdz_int,1),std(dEdz_int),'-b','linewidth',3.0);
% plot(Efit,dEdz_2,'-k','linewidth',6.0);
% plot(scattdata.optical.inel_dcsdata.E,scattdata.optical.inel_dcsdata.Sp,'-r','linewidth',3.0);
plot(scattdata.optical.E,scattdata.optical.Sp,'-r','linewidth',3.0);
% if exist('TruongData')
%     plot(TruongData.PMMA.E,TruongData.PMMA.Sp,'--ob');
% end

xlabel('E (eV)','fontsize',30);
ylabel('dE/dz (eV/nm)','fontsize',30);
set(gca,'XScale','log','YScale','linear','fontsize',30,'linewidth',3.0);

%% Acid image post-processing
close all;

pagsum_zx=sum(pagimg_orig,3);
pagsum_yx=reshape(sum(pagimg_orig,1),[size(pagimg_orig,3) size(pagimg_orig,2)]);
pagsum_zy=reshape(sum(pagimg_orig,2),[size(pagimg_orig,1) size(pagimg_orig,3)]);

acidsum_zx=sum(acidimg,3);
acidsum_yx=reshape(sum(acidimg,1),[size(acidimg,3) size(acidimg,2)]);
acidsum_zy=reshape(sum(acidimg,2),[size(acidimg,1) size(acidimg,3)]);

psf=fspecial('gaussian',[100 100],7.5/0.2);
% acidsum_zx=conv2(acidsum_zx,psf,'same');
% acidsum_yx=conv2(acidsum_yx,psf,'same');
% acidsum_zy=conv2(acidsum_zy,psf,'same');

zdist=acid_xyz(:,3)-zval; % distance in z;
zdist2=zdist;
% zdist2(zdist2<=zval)=0;
zdist2=zdist2(zdist2>0);
depth_mean=mean(zdist2);
depth_std=std(zdist2);

fprintf('Depth Mean = %.4f nm;',depth_mean);
fprintf('  Depth Stdev = %.4f nm\n',depth_std);

nAcids_thruz=[];
zvec=[];
for zcount=1:size(acidimg,1) % go through each slice in z-coordinate, get the number of acids
    acidimg_tmp=acidimg(zcount,:,:);
    nAcids_thruz(zcount)=sum(acidimg_tmp(:));
    zvec(zcount)=univ.grid.z(zcount,1,1);
end

% nAcids_thruz(zvec<=zval)=0;

nAcids_thruz=nAcids_thruz(zvec>=zval);
zvec=zvec(zvec>=zval);
zvec=zvec-zval;

nAcids_norm=nAcids_thruz./sum(nAcids_thruz);
nAcids_cumulsum=cumsum(nAcids_norm);

fprintf('# of incident electrons = %d; # of Acids created = %.4f; QE = %.4f\n',sum(elecimg_inc(:)),sum(nAcids_thruz),sum(nAcids_thruz)/sum(elecimg_inc(:)));
% depth_mean=sum(zvec.*nAcids_norm);
% depth_std=sum((zvec-depth_mean).^2.*nAcids_norm);

% figure;
hold on;
plot(zvec,nAcids_thruz,'-ok');
xlabel('z (nm)');
ylabel('# of Acids');

figure;plot(zvec,nAcids_cumulsum,'-ob');
xlabel('z (nm)');
ylabel('Cumulative Sum');

figure;imagesc(acidsum_zx);title('Acid zx');
figure;imagesc(acidsum_yx);title('Acid yx');
figure;imagesc(pagsum_zx);title('Pre-Exose PAG zx');
figure;imagesc(pagsum_yx);title('Pre-Expose PAG yx');

%% save the workspace for use in the low-energy contrast curve simulator
clc;
close all;

save(sprintf('LEEM_SimResults_2\\Ein=%.2f_Dose=%.2f_Ef=15.5_pag-Emin=5_rcnrad=2_PAG=0.4_T2.mat',Esweep,Dose));
% save(sprintf('LEEM_SimResults_2\\Ein=%.2f_Dose=%.2f_Ef=15.5_pag-Emin=5_rcnrad=2_PAG=0.4_20nmX20nm_T1.mat',Esweep,Dose));

% save(sprintf('LEEM_SimResults\\Ein=%.2f_Dose=%.2f_Ef=15.5_Emin=3_rcnrad=2_PAG=0.4_24nmArea_PhononModel=Frohlich.mat',Esweep,Dose));

%% LEEM thru-dose analysis
clear;
clc;
close all;

% data20=load('LEEM_SimResults\Ein=20eV_10X10EBeam_WithVibr.mat');
% data40=load('LEEM_SimResults\Ein=40eV_10X10EBeam_WithVibr.mat');
% data40_d0p1=load('LEEM_SimResults\Ein=40eV_Dose=0p1_10X10EBeam_WithVibr.mat');
% data40_d0p5=load('LEEM_SimResults\Ein=40eV_Dose=0p5_10X10EBeam_WithVibr.mat');
% data40_d1=load('LEEM_SimResults\Ein=40eV_Dose=1_10X10EBeam_WithVibr.mat');
% data40_d1p5=load('LEEM_SimResults\Ein=40eV_Dose=1p5_10X10EBeam_WithVibr.mat');
% data40_d2=load('LEEM_SimResults\Ein=40eV_Dose=2_10X10EBeam_WithVibr.mat');

% data80_d0p1=load('LEEM_SimResults\Ein=80eV_Dose=0p1_10X10EBeam_WithVibr.mat');
% data80_d0p5=load('LEEM_SimResults\Ein=80eV_Dose=0p5_10X10EBeam_WithVibr.mat');
% data80_d1=load('LEEM_SimResults\Ein=80eV_Dose=1_10X10EBeam_WithVibr.mat');
% data80_d1p5=load('LEEM_SimResults\Ein=80eV_Dose=1p5_10X10EBeam_WithVibr.mat');
% data80_d2=load('LEEM_SimResults\Ein=80eV_Dose=2_10X10EBeam_WithVibr.mat');
% data80=load('LEEM_SimResults\Ein=80eV_10X10nm-EBeam_WithVibr.mat');

% data{1}=load('LEEM_SimResults\Ein=80eV_Dose=1_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=0.5nm.mat');
% data{2}=load('LEEM_SimResults\Ein=80eV_Dose=1_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=0.5nm.mat');
% data{3}=load('LEEM_SimResults\Ein=80eV_Dose=1_RCNRAD=3nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=0.5nm.mat');
% data{4}=load('LEEM_SimResults\Ein=80eV_Dose=1_RCNRAD=5nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=0.5nm.mat');
% data{5}=load('LEEM_SimResults\Ein=80eV_Dose=5_RCNRAD=5nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=0.5nm.mat');

% data{1}=load('LEEM_SimResults\Ein=80eV_Dose=0p1_10X10EBeam_WithVibr.mat');
% data{2}=load('LEEM_SimResults\Ein=80eV_Dose=0p5_10X10EBeam_WithVibr.mat');
% data{3}=load('LEEM_SimResults\Ein=80eV_Dose=1_10X10EBeam_WithVibr.mat');
% data{4}=load('LEEM_SimResults\Ein=80eV_Dose=1p5_10X10EBeam_WithVibr.mat');
% data{5}=load('LEEM_SimResults\Ein=80eV_Dose=2_10X10EBeam_WithVibr.mat');

data={};
data{length(data)+1}=load('LEEM_SimResults\Ein=40eV_Dose=0.5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=40eV_Dose=1_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=40eV_Dose=2_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=40eV_Dose=3_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=40eV_Dose=5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=40eV_Dose=10_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=40eV_Dose=20_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=80eV_Dose=1_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=80eV_Dose=2_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=80eV_Dose=3_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=80eV_Dose=5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=80eV_Dose=10_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
data{length(data)+1}=load('LEEM_SimResults\Ein=80eV_Dose=20_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');

%% post-proc of LEEM analysis results: RCN/Diffn/Base Model
clc;
close all;

base_load=0.1; % PAG is 0.5/nm3;
acid_difflen=10; % nm, over a 90 s PEB time at 100 deg. C.
base_difflen=10; % nm, over a 90 s PEB time at 100 deg. C.
prot_load=1; % /nm3
kD=5; % Deprotection Rate (nm3/s)
kQ=15; % A/B Quenching Rate (nm3/s)
tPEB=90;
dt_PEB=1;
time=dt_PEB:dt_PEB:tPEB;

for i = 1:length(data)
% for i = 1
    fprintf('Sweep Data %d of %d\n',i,length(data));
    acidimg=data{i}.acidimg;
    rng(102);
    baseimg=base_load.*ones(size(acidimg)).*prod(data{i}.univ.px_nm);
    baseimg=poissrnd(baseimg);
    protimg=prot_load.*ones(size(acidimg)).*prod(data{i}.univ.px_nm);
    protimg=poissrnd(protimg);
    protimg_pre=protimg;
    rng('shuffle');
    
    %%%% define the 3D PSF:
    xgrid=data{i}.univ.grid.x;
    ygrid=data{i}.univ.grid.y;
    zgrid=data{i}.univ.grid.z;
    Acidsig=acid_difflen/sqrt(length(time));
    Basesig=base_difflen/sqrt(length(time));
    
    psf_radius=sqrt(xgrid.^2+ygrid.^2+zgrid.^2);
    acidpsf=1/(Acidsig^3*(2*pi)^(3/2)).*exp(-psf_radius.^2./(2*Acidsig^2));
    basepsf=1/(Basesig^3*(2*pi)^(3/2)).*exp(-psf_radius.^2./(2*Basesig^2));
    
    %%%% convert rcn rates to pixels/s
    kD_px=kD*prod(data{i}.univ.px_nm);
    kQ_px=kQ*prod(data{i}.univ.px_nm);
    
    
    for j = 1:length(time)
        fprintf('...PEB Time %d of %d\n',j,length(time));
        tstart_PEBStep=tic;
        acidimg=convn(acidimg,acidpsf,'same');
        baseimg=convn(baseimg,basepsf,'same');
        Nrcn_AB=kQ_px.*acidimg.*baseimg;
        Ndeprot=kD_px.*acidimg.*protimg;
        
        protimg=protimg-Ndeprot;protimg(protimg<0)=0;
        acidimg=acidimg-Nrcn_AB;acidimg(acidimg<0)=0;
        baseimg=baseimg-Nrcn_AB;baseimg(baseimg<0)=0;
        tend_PEBStep=toc(tstart_PEBStep);
        fprintf('......Took %.4f s\n',tend_PEBStep);
    end
    deprimg=(protimg_pre-protimg)./protimg_pre;
    deprimg(isnan(deprimg))=0;
    deprimg(deprimg==Inf)=0;
    deprimg(deprimg==-Inf)=0;
    
end

%% post-proc of LEEM analysis results [generalized sweep]: 1-D convolution
close all;
clc;
psfsig=5/2.35;

fig1=figure;
hold on;

fig2=figure;
hold on;

% fig3=figure;
% hold on;

% fig_styles={'-ob','-or','-ok','-og','--ob','--or','--ok','--og'};
fig_styles={'-b','-r','-k','-g','--b','--r','--k','--g'};

fig_styles=repmat(fig_styles,[1,100]);

CD_thr=2; % the absolute threshold for CD determination
CD=[];

for i = 1:length(data)
    z1=data{i}.zvec;
    z=[sort(-z1) z1];
    
    psf=1/(psfsig*sqrt(2*pi)).*exp(-(z).^2./(2*psfsig^2));
    
    acids1=data{i}.nAcids_thruz;
    acids=interp1(z1,data{i}.nAcids_thruz,z);
    acids(isnan(acids))=0;
    depr=conv(acids,psf,'same');
    depr_max=max(depr);
    numacids(i)=sum(data{i}.acidimg(:));
    
    z2=z(z>=0);
    depr2=depr(z>=0);
    
    if CD_thr>=min(depr2) & CD_thr<=max(depr2)
        idx=find(depr2==CD_thr);
        if isempty(idx)
            idx1=find(depr2<CD_thr);idx1=idx1(1);
            idx2=find(depr2>CD_thr);idx2=idx2(end);
            p=polyfit([z2(idx1) z2(idx2)],[depr2(idx1) depr2(idx2)],1);
            CD(i)=(CD_thr-p(2))/p(1);
        else
            CD(i)=z2(idx);
        end
    else
        CD(i)=0;
    end
    
    
%     figure(fig3);
%     plot(z,z,fig_styles{i});
    
    figure(fig1);
    plot(z,acids./mean(acids1(1:5)),fig_styles{i});
    xlim([0 max(z)]);
    
    figure(fig2);
%     plot(z,depr./max(depr),fig_styles{i},'linewidth',3.0);
    plot(z,depr,fig_styles{i},'linewidth',3.0);
    xlim([0 max(z)]);
    grid on;
end

%% post-proc of LEEM analysis results [when swept parameter is dose]

dosevec=[0.5 1 2 3 5 10 20];

figure;
plot(dosevec,50-CD(1:length(dosevec)),'-ob');
hold on;
plot(dosevec,50-CD(length(dosevec)+1:length(CD)),'-or');
set(gca,'XScale','linear');

%% post-proc of LEEM thru-dose analysis [OLD, very specific]
close all;
psfsig=15/2.35;

% psf20=1/(psfsig*sqrt(2*pi)).*exp(-(data20.zvec).^2./(2*psfsig^2));

z40=data40_d0p5.zvec;
z40=[sort(-z40) z40];
psf40=1/(psfsig*sqrt(2*pi)).*exp(-(z40).^2./(2*psfsig^2));

% psf80=1/(psfsig*sqrt(2*pi)).*exp(-(data80.zvec).^2./(2*psfsig^2));

% depr20=conv(data20.nAcids_thruz,psf20,'same');

%%%%% 80 eV:
doses_E80=[0.1 0.5 1 1.5 2];
acids=interp1(data80_d0p1.zvec,data80_d0p1.nAcids_thruz,z40);
acids(isnan(acids))=0;
depr80_d0p1=conv(acids,psf40,'same');
depr80_d0p1_max=max(depr80_d0p1);
numacids_d80(1)=sum(data80_d0p1.acidimg(:));

acids=interp1(data80_d0p5.zvec,data80_d0p5.nAcids_thruz,z40);
acids(isnan(acids))=0;
depr80_d0p5=conv(acids,psf40,'same');
depr80_d0p5_max=max(depr80_d0p5);
numacids_d80(2)=sum(data80_d0p5.acidimg(:));

acids=interp1(data80_d1.zvec,data80_d1.nAcids_thruz,z40);
acids(isnan(acids))=0;
depr80_d1=conv(acids,psf40,'same');
depr80_d1_max=max(depr80_d1);
numacids_d80(3)=sum(data80_d1.acidimg(:));

acids=interp1(data80_d1p5.zvec,data80_d1p5.nAcids_thruz,z40);
acids(isnan(acids))=0;
depr80_d1p5=conv(acids,psf40,'same');
depr80_d1p5_max=max(depr80_d1p5);
numacids_d80(4)=sum(data80_d1p5.acidimg(:));

acids=interp1(data80_d2.zvec,data80_d2.nAcids_thruz,z40);
acids(isnan(acids))=0;
depr80_d2=conv(acids,psf40,'same');
depr80_d2_max=max(depr80_d2);
numacids_d80(5)=sum(data80_d2.acidimg(:));

%%%%% 40 eV:
doses_E40=[0.1 0.5 1 1.5 2];
acids=interp1(data40_d0p1.zvec,data40_d0p1.nAcids_thruz,z40);
acids(isnan(acids))=0;
depr40_d0p1=conv(acids,psf40,'same');
depr40_d0p1_max=max(depr40_d0p1);
numacids_d40(1)=sum(data40_d0p1.acidimg(:));

acids=interp1(data40_d0p5.zvec,data40_d0p5.nAcids_thruz,z40);
acids(isnan(acids))=0;
depr40_d0p5=conv(acids,psf40,'same');
depr40_d0p5_max=max(depr40_d0p5);
numacids_d40(2)=sum(data40_d0p5.acidimg(:));

acids=interp1(data40_d1.zvec,data40_d1.nAcids_thruz,z40);
acids(isnan(acids))=0;
depr40_d1=conv(acids,psf40,'same');
depr40_d1_max=max(depr40_d1);
numacids_d40(3)=sum(data40_d1.acidimg(:));

acids=interp1(data40_d1p5.zvec,data40_d1p5.nAcids_thruz,z40);
acids(isnan(acids))=0;
depr40_d1p5=conv(acids,psf40,'same');
depr40_d1p5_max=max(depr40_d1p5);
numacids_d40(4)=sum(data40_d1p5.acidimg(:));

acids=interp1(data40_d2.zvec,data40_d2.nAcids_thruz,z40);
acids(isnan(acids))=0;
depr40_d2=conv(acids,psf40,'same');
depr40_d2_max=max(depr40_d2);
numacids_d40(5)=sum(data40_d2.acidimg(:));

% depr80=conv(data80.nAcids_thruz,psf80,'same');

figure;
% plot(data20.zvec,data20.nAcids_thruz,'-b','linewidth',3.0);
hold on; 
plot(data40_d0p1.zvec,data40_d0p1.nAcids_thruz,'-b','linewidth',3.0);
plot(data40_d0p5.zvec,data40_d0p5.nAcids_thruz,'-r','linewidth',3.0);
plot(data40_d1.zvec,data40_d1.nAcids_thruz,'-k','linewidth',3.0);
plot(data40_d1p5.zvec,data40_d1p5.nAcids_thruz,'-g','linewidth',3.0);
plot(data40_d2.zvec,data40_d2.nAcids_thruz,'--y','linewidth',3.0);

plot(data80_d0p1.zvec,data80_d0p1.nAcids_thruz,'--b','linewidth',3.0);
plot(data80_d0p5.zvec,data80_d0p5.nAcids_thruz,'--r','linewidth',3.0);
plot(data80_d1.zvec,data80_d1.nAcids_thruz,'--k','linewidth',3.0);
plot(data80_d1p5.zvec,data80_d1p5.nAcids_thruz,'--g','linewidth',3.0);
plot(data80_d2.zvec,data80_d2.nAcids_thruz,'--y','linewidth',3.0);

% plot(data80.zvec,data80.nAcids_thruz,'-k','linewidth',3.0);
xlabel('z (nm)','fontsize',30);
ylabel('# of Acids','fontsize',30);
set(gca,'fontsize',30);

% figure;
% % plot(data20.zvec,data20.nAcids_thruz,'-b','linewidth',3.0);
% hold on; 
% plot(data40_d0p1.zvec,data40_d0p1.nAcids_thruz./max(data40_d0p1.nAcids_thruz),'-b','linewidth',3.0);
% plot(data40_d0p5.zvec,data40_d0p5.nAcids_thruz./max(data40_d0p5.nAcids_thruz),'-r','linewidth',3.0);
% plot(data40_d1.zvec,data40_d1.nAcids_thruz./max(data40_d1.nAcids_thruz),'-k','linewidth',3.0);
% plot(data40_d1p5.zvec,data40_d1p5.nAcids_thruz./max(data40_d1p5.nAcids_thruz),'-g','linewidth',3.0);
% % plot(data80.zvec,data80.nAcids_thruz,'-k','linewidth',3.0);
% xlabel('z (nm)','fontsize',30);
% ylabel('Deprotection Profile','fontsize',20);

figure;
% plot(data20.zvec,depr20,'b','linewidth',3.0);
hold on;
plot(z40,depr40_d0p1./max(depr40_d0p1).*max(depr80_d0p1),'b','linewidth',3.0);
plot(z40,depr40_d0p5./max(depr40_d0p5).*max(depr80_d0p5),'r','linewidth',3.0);
plot(z40,depr40_d1./max(depr40_d1).*max(depr80_d1),'k','linewidth',3.0);
plot(z40,depr40_d1p5./max(depr40_d1p5).*max(depr80_d1p5),'g','linewidth',3.0);
plot(z40,depr40_d2./max(depr40_d2).*max(depr80_d2),'-y','linewidth',3.0);

plot(z40,depr80_d0p1,'--b','linewidth',3.0);
plot(z40,depr80_d0p5,'--r','linewidth',3.0);
plot(z40,depr80_d1,'--k','linewidth',3.0);
plot(z40,depr80_d1p5,'--g','linewidth',3.0);
plot(z40,depr80_d2,'--y','linewidth',3.0);

xlim([0 30]);
xlabel('z (nm)','fontsize',30);
ylabel('Deprotection Profile','fontsize',30);
set(gca,'fontsize',30);

figure;
% plot(data20.zvec,depr20,'b','linewidth',3.0);
hold on;
plot(z40,depr40_d0p1,'b','linewidth',3.0);
plot(z40,depr40_d0p5,'r','linewidth',3.0);
plot(z40,depr40_d1,'k','linewidth',3.0);
plot(z40,depr40_d1p5,'g','linewidth',3.0);
plot(z40,depr40_d2,'-y','linewidth',3.0);

plot(z40,depr80_d0p1,'--b','linewidth',3.0);
plot(z40,depr80_d0p5,'--r','linewidth',3.0);
plot(z40,depr80_d1,'--k','linewidth',3.0);
plot(z40,depr80_d1p5,'--g','linewidth',3.0);
plot(z40,depr80_d2,'--y','linewidth',3.0);

xlim([0 30]);
xlabel('z (nm)','fontsize',30);
ylabel('Deprotection Profile','fontsize',30);
set(gca,'fontsize',30);
legend('Dose=0.1 e^-/pixel','Dose=0.5 e^-/pixel','Dose=1 e^-/pixel','Dose=1.5 e^-/pixel','Dose=2 e^-/pixel');

figure;
plot(doses_E40,[depr40_d0p1_max depr40_d0p5_max depr40_d1_max depr40_d1p5_max depr40_d2_max],'-ob','linewidth',3.0);
hold on;
plot(doses_E80,[depr80_d0p1_max depr80_d0p5_max depr80_d1_max depr80_d1p5_max depr80_d2_max],'-or','linewidth',3.0);
xlabel('Dose (e^-/pixel)','fontsize',30);
ylabel('Max. Deprot. Value','fontsize',30);
set(gca,'fontsize',30);

%% Thru-trial acid stats: Get distribution of mean radii
clc;
close all;

radius_stats=[];
radcount=1;
for i = 1:length(radius_acids_thrutrial)
    if ~isnan(radius_acids_thrutrial{i})
        radvec(radcount)=mean(radius_acids_thrutrial{i});
        radcount=radcount+1;
    end
end

fprintf('<radius> = %.4f; sigma_radius = %.4f\n',mean(radvec),std(radvec));
fprintf('<Acids> = %.4f; sigma_Acids = %.4f\n\n',mean(nacids_thrutrial),std(nacids_thrutrial));

figure;hist(radvec,100);
xlabel('Radius (nm)');
ylabel('Count');

%% Mean Free Path Components Analysis
clc;
close all;

na=1.2/120*6*1e23;

%%%% vibrational data
vibr_E=scattdata.vibr.ics(:,1);
vibr_ics=scattdata.vibr.ics(:,2:end);
vibr_mfp=1./(vibr_ics.*scattdata.vibr.ics_mult.*na).*1e7;

%%%% optical data
opt_E=scattdata.optical.E;
opt_imfp=scattdata.optical.imfp;

%%%% Fit
Efit_min=min([reshape(vibr_E,[1,length(vibr_E)]) reshape(opt_E,[1,length(opt_E)])]);
Efit_max=max([reshape(vibr_E,[1,length(vibr_E)]) reshape(opt_E,[1,length(opt_E)])]);
Efit=Efit_min:Efit_max;

mfp_fit=[];
for i = 1:length(Efit)
    for j = 1:size(vibr_mfp,2)
        tmpvec=vibr_mfp(:,j);
        idx=find(~isnan(tmpvec)==1);
        mfp_fit(i,j)=interp1(vibr_E(idx),tmpvec(idx),Efit(i));
    end
    mfp_fit(i,j+1)=interp1(opt_E,opt_imfp,Efit(i));
end
mfp_fit(isnan(mfp_fit))=Inf;
mfp_fit_recip=1./mfp_fit;
mfp_fit2=1./sum(mfp_fit_recip,2);

figure;
plot(vibr_E,vibr_mfp(:,1),'-b','linewidth',6.0);
hold on;
plot(vibr_E,vibr_mfp(:,2),'-r','linewidth',6.0);
plot(vibr_E,vibr_mfp(:,3),'-k','linewidth',6.0);
plot(vibr_E,vibr_mfp(:,4),'-g','linewidth',6.0);
plot(opt_E,opt_imfp,'--b','linewidth',6.0);
plot(Efit',mfp_fit2,'--r','linewidth',6.0);
xlabel('Energy (eV)','fontsize',30);
ylabel('Mean Free Path (nm)','fontsize',30);
set(gca,'YScale','log','XScale','log','fontsize',30,'linewidth',3.0);
xlim([0 200]);
legend('0.173 eV Mode','0.353 eV Mode','0.531 eV Mode','0.7 eV Mode','Inel. MFP','Total MFP');

figure;
plot(scattdata.vibr.ics(:,1),scattdata.vibr.ics(:,2).*scattdata.vibr.ics_mult,'-ob','linewidth',3.0,'MarkerSize',12)
hold on;
plot(scattdata.vibr.ics(:,1),scattdata.vibr.ics(:,3).*scattdata.vibr.ics_mult,'-or','linewidth',3.0,'MarkerSize',12)
plot(scattdata.vibr.ics(:,1),scattdata.vibr.ics(:,4).*scattdata.vibr.ics_mult,'-ok','linewidth',3.0,'MarkerSize',12)
plot(scattdata.vibr.ics(:,1),scattdata.vibr.ics(:,5).*scattdata.vibr.ics_mult,'-og','linewidth',3.0,'MarkerSize',12)
set(gca,'XScale','linear','YScale','linear','fontsize',30);
xlabel('Incident Electron Energy (eV)','fontsize',30);
ylabel('TCS (cm^2)','fontsize',30);
legend('0.173 eV Mode','0.353 eV Mode','0.531 eV Mode','0.7 eV Mode');

xlim_vec=[0 130];
ylim_vec=[0 0.35];
fig_fontsize=30;
fig_lw=3.0;

figure;
subplot(3,2,1);
plot(scattdata.vibr.Epr1.angledata(:,1),scattdata.vibr.Epr1.angledata(:,2),'b','linewidth',6.0);
hold on;
plot(scattdata.vibr.Epr1.angledata(:,1),scattdata.vibr.Epr1.angledata(:,3),'r','linewidth',6.0);
% title(sprintf('Eo = %.2f eV',scattdata.vibr.Epr1.Eo));
xlim(xlim_vec);
% ylim(ylim_vec);
set(gca,'fontsize',fig_fontsize,'linewidth',fig_lw);

subplot(3,2,2);
plot(scattdata.vibr.Epr2.angledata(:,1),scattdata.vibr.Epr2.angledata(:,2),'b','linewidth',6.0);
hold on;
plot(scattdata.vibr.Epr2.angledata(:,1),scattdata.vibr.Epr2.angledata(:,3),'r','linewidth',6.0);
plot(scattdata.vibr.Epr2.angledata(:,1),scattdata.vibr.Epr2.angledata(:,4),'k','linewidth',6.0);
% title(sprintf('Eo = %.2f eV',scattdata.vibr.Epr2.Eo));
xlim(xlim_vec);
% ylim(ylim_vec);
set(gca,'fontsize',fig_fontsize,'linewidth',fig_lw);

subplot(3,2,3);
plot(scattdata.vibr.Epr3.angledata(:,1),scattdata.vibr.Epr3.angledata(:,2),'b','linewidth',6.0);
hold on;
plot(scattdata.vibr.Epr3.angledata(:,1),scattdata.vibr.Epr3.angledata(:,3),'r','linewidth',6.0);
plot(scattdata.vibr.Epr3.angledata(:,1),scattdata.vibr.Epr3.angledata(:,4),'k','linewidth',6.0);
plot(scattdata.vibr.Epr3.angledata(:,1),scattdata.vibr.Epr3.angledata(:,5),'g','linewidth',6.0);
% title(sprintf('Eo = %.2f eV',scattdata.vibr.Epr3.Eo));
xlim(xlim_vec);
% ylim(ylim_vec);
set(gca,'fontsize',fig_fontsize,'linewidth',fig_lw);

subplot(3,2,4);
plot(scattdata.vibr.Epr4.angledata(:,1),scattdata.vibr.Epr4.angledata(:,2),'b','linewidth',6.0);
hold on;
plot(scattdata.vibr.Epr4.angledata(:,1),scattdata.vibr.Epr4.angledata(:,3),'r','linewidth',6.0);
plot(scattdata.vibr.Epr4.angledata(:,1),scattdata.vibr.Epr4.angledata(:,4),'k','linewidth',6.0);
plot(scattdata.vibr.Epr4.angledata(:,1),scattdata.vibr.Epr4.angledata(:,5),'g','linewidth',6.0);
% title(sprintf('Eo = %.2f eV',scattdata.vibr.Epr4.Eo));
xlim(xlim_vec);
% ylim(ylim_vec);
set(gca,'fontsize',fig_fontsize,'linewidth',fig_lw);

subplot(3,2,5);
plot(scattdata.vibr.Epr5.angledata(:,1),scattdata.vibr.Epr5.angledata(:,2),'b','linewidth',6.0);
hold on;
plot(scattdata.vibr.Epr5.angledata(:,1),scattdata.vibr.Epr5.angledata(:,3),'r','linewidth',6.0);
plot(scattdata.vibr.Epr5.angledata(:,1),scattdata.vibr.Epr5.angledata(:,4),'k','linewidth',6.0);
plot(scattdata.vibr.Epr5.angledata(:,1),scattdata.vibr.Epr5.angledata(:,5),'g','linewidth',6.0);
% title(sprintf('Eo = %.2f eV',scattdata.vibr.Epr5.Eo));
xlim(xlim_vec);
% ylim(ylim_vec);
set(gca,'fontsize',fig_fontsize,'linewidth',fig_lw);

subplot(3,2,6);
plot(scattdata.vibr.Epr6.angledata(:,1),scattdata.vibr.Epr6.angledata(:,2),'b','linewidth',6.0);
hold on;
plot(scattdata.vibr.Epr6.angledata(:,1),scattdata.vibr.Epr6.angledata(:,3),'r','linewidth',6.0);
plot(scattdata.vibr.Epr6.angledata(:,1),scattdata.vibr.Epr6.angledata(:,4),'k','linewidth',6.0);
plot(scattdata.vibr.Epr6.angledata(:,1),scattdata.vibr.Epr6.angledata(:,5),'g','linewidth',6.0);
% title(sprintf('Eo = %.2f eV',scattdata.vibr.Epr6.Eo));
xlim(xlim_vec);
% ylim(ylim_vec);
set(gca,'fontsize',fig_fontsize,'linewidth',fig_lw);

%% Mean Free Path Components Analysis
clc;
% close all;

Eopt=scattdata.optical.E;
imfp_opt=scattdata.optical.imfp;

E_interp=1:1:max(Eopt);

imfp_opt_interp=interp1(Eopt,imfp_opt,E_interp);
imfp_opt_interp(E_interp<min(Eopt))=Inf;

% E,eps0,epsinf,E_optphonon
imfp_vibr=scattdata.vibr.imfp_func(E_interp,scattdata.vibr.eps0,scattdata.vibr.epsinf,scattdata.vibr.hbarw);
% imfp_vibr=scattdata.vibr.imfp_func(E_interp,scattdata.vibr.eps0,scattdata.vibr.epsinf,0.2);
imfp_vibr1=scattdata.vibr.imfp_func(E_interp,scattdata.vibr.eps0,scattdata.vibr.epsinf,0.102);
imfp_vibr2=scattdata.vibr.imfp_func(E_interp,scattdata.vibr.eps0,scattdata.vibr.epsinf,0.109);
imfp_vibr=1./(1./imfp_vibr1 + 1./imfp_vibr2);

imfp_tot=1./(1./imfp_vibr + 1./imfp_opt_interp);

figure;
plot(E_interp,imfp_opt_interp,'b','linewidth',3.0);
hold on;
plot(E_interp,imfp_vibr,'k','linewidth',3.0);
plot(E_interp,imfp_tot,'r','linewidth',3.0);
ylim([0 28]);
xlabel('Electron Energy (eV)');
ylabel('Inelastic Mean Free Path (nm)');
set(gca,'fontsize',30,'linewidt',3.0);
% legend('Dielectric Model','Electron-Phonon Frohlich Model','Combined');

%% ANgle and Energy distributions: Dielectric Model
clc;
close all;

Esel_idx=[17 27 37];

plt_styles={'-b','-k','-r','--b','--k','--r'};

%%%% ANgle Distribution
legstr={};
imfp_theta=[];
theta_mean=[];
figure;hold on;
for i = 1:length(Esel_idx)
    xvec=scattdata.optical.inel_dcsdata.thetamat(Esel_idx(i),:);
    yvec=scattdata.optical.inel_dcsdata.dsigdOmega(Esel_idx(i),:);
    imfp_theta(i)=1/trapz(xvec,2*pi.*sin(xvec).*yvec.*1e-9);
    theta_mean(i)=trapz(xvec,xvec.*yvec)/trapz(xvec,yvec);
    plot(xvec.*180/pi,yvec.*1e-9,plt_styles{i},'linewidth',3.0);
    legstr{i}=sprintf('E = %.2f eV',scattdata.optical.inel_dcsdata.E(Esel_idx(i)));
end
legend(legstr);
xlabel('Theta (Degrees)');
ylabel('d\sigma/d\Omega');
xlim([0 90]);
set(gca,'linewidth',3.0,'fontsize',30);
box on;

%%%% ENergy Distribution
legstr={};
imfp_E=[];
Eloss_mean=[];
figure;hold on;
for i = 1:length(Esel_idx)
    xvec=scattdata.optical.inel_dcsdata.Elossmat(Esel_idx(i),:);
    yvec=scattdata.optical.inel_dcsdata.dsigdE(Esel_idx(i),:);
    imfp_E(i)=1/trapz(xvec,yvec.*1e-9);
    Eloss_mean(i)=trapz(xvec,xvec.*yvec)/trapz(xvec,yvec);
    plot(xvec,yvec.*1e-9,plt_styles{i},'linewidth',3.0);
    legstr{i}=sprintf('E = %.2f eV',scattdata.optical.inel_dcsdata.E(Esel_idx(i)));
end
legend(legstr);
xlabel('Eloss (eV)');
ylabel('d\sigma/d\E_{loss}');
set(gca,'linewidth',3.0,'fontsize',30);
box on;
%% angle scattering distribution of phonon-electron process
clc;
close all;

ntrials=10000;
El_Ph_Angle_Dist(scattdata.vibr,ntrials,10);


%% Scattering Sim: Optical + Vibrational Excitations OLD (Likely deletable) [Real Resist 3D Grid]

clear;
clc;
close all;

addpath('F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\GlobalFunctions');

scattdata.vibr=load('VibrExcit_Data_Khakoo_2013.mat');
% scattdata.vibr=load('Phenol_Vibr+Electronic_Excit_Data_FlindersUniv_2015.mat');

optdata_path='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\IMFP_Analysis\MatFiles\IMFP_Components\';
scattdata.optical=load([optdata_path 'Sp_IMFP_Inelastic_Components_Ef=5eV_Eloss=[0 500]_v2.mat']); % _v2 has data structures better suited to the subsequently called programs

inel_dcsdata_path='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\Traj_MonteCarlo\';
scattdata.optical.inel_dcsdata=load([inel_dcsdata_path 'Inelastic_DCS_Data.mat']);

%%%%%% create PAG grid
rho_pag=0.5; % per nm3
px_nm=[0.1 0.1 0.1];

univ.size_nm=[10 10 10]; % nm
univ.px_nm=px_nm;
univ.npx=round(univ.size_nm./univ.px_nm);

xgrid=-univ.size_nm(1)/2:px_nm(1):univ.size_nm(1)/2-px_nm(1);
ygrid=-univ.size_nm(2)/2:px_nm(2):univ.size_nm(2)/2-px_nm(2);
zgrid=-univ.size_nm(3)/2:px_nm(3):univ.size_nm(3)/2-px_nm(3);
[xgrid,ygrid,zgrid]=meshgrid(xgrid,ygrid,zgrid);
univ.grid.x=xgrid;
univ.grid.y=ygrid;
univ.grid.z=zgrid;

pagimg=rho_pag*prod(univ.px_nm).*ones(univ.npx);
pagimg=poissrnd(pagimg);
elecimg=zeros(size(pagimg));

%%%% Define the exposure mask
mask=zeros(univ.npx(1:2));
mask(40:60,40:60,:)=1;

%%%% Generate LP-FIltered image
psf=fspecial('gaussian',[10,10],5);
aimg=imfilter(mask,psf);
% figure;imagesc(aimg);colormap('gray');colorbar;title('aerial image');

%%%% generate 3-D photon absorption image
Dose=50; % mJ/cm2
Ephoton=6.626*1e-34*3*1e8/(13.5*1e-9);
pixel_area=prod(univ.px_nm(1:2))*1e-14; % convert to cm2
nph_perpixel=Dose*0.001/Ephoton*pixel_area;

% pimg=aimg.*nph_perpixel;
% pimg=poissrnd(pimg);
% figure;imagesc(pimg);colormap('gray');colorbar;title('incident photon image');

absorptivity=0.0042;
% absimg=zeros([100 100 100]);
% absimg=photon_abs(pimg,absorptivity,absimg);
% figure;imagesc(sum(absimg,3));colormap('gray');colorbar;title('absorbed photon image');

%%%%%% Initialize the event structure [holds info about excitations triggered]
% event{1}.xyz=round(univ.npx/2);
event{1}.xyz=[0 0 0];
xyzglobal.x=event{1}.xyz(1);
xyzglobal.y=event{1}.xyz(2);
xyzglobal.z=event{1}.xyz(3);

event{1}.Ein=80;
event{1}.Eout=0;
event{1}.Eloss=10;
event{1}.Ese=80;
event{1}.act='SE';
event{1}.nacid=0;
event{1}.nSE=0;
event{1}.elecimg=elecimg;
event{1}.pag.img=pagimg;
event{1}.pag.rcnrad=2; % nm
event{1}.pag.rho=rho_pag;
event{1}.controlparms.model_pag_devel=1;
event{1}.univ=univ;

% ntrials=sum(absimg(:));
% abspos=find(absimg~=0);
ntrials=200;

plotfig=0;
if plotfig~=0
    fig1=figure;hold on;colors='brgk';
end

radius_acids=[];
nacids_thrutrial=[];
tstart=tic;

for trial_count=1:ntrials
    tst_trials=tic;
    fprintf('Trial %d of %d\n',trial_count,ntrials);
%     nphotons=absimg(abspos(trial_count));
    nphotons=1;
%     [row,col,slice]=ind2sub(size(absimg),abspos(trial_count));
%     x_abs=univ.grid.x(row,col,slice);
%     y_abs=univ.grid.y(row,col,slice);
%     z_abs=univ.grid.z(row,col,slice);
    for photon_count=1:nphotons
        event{1}.xyz=[0 0 0];
        xyzglobal.x=event{1}.xyz(1);
        xyzglobal.y=event{1}.xyz(2);
        xyzglobal.z=event{1}.xyz(3);

        eventdata={};
        [eventdata]=Scattcalc_lowE(event,scattdata,eventdata,xyzglobal);

        Elossvec=[];act={};nSE=[];imfp=[];theta=[];phi=[];xyz=[];nacid=[];
        xyz_acids=[];%radius_acids=[];
        for count = 1:length(eventdata)
            Elossvec(count)=eventdata{count}.Eloss;
            xyz(count,:)=eventdata{count}.xyz;
            act{count}=eventdata{count}.act;
            nSE(count)=eventdata{count}.nSE;
            nacids(count)=eventdata{count}.nacid;
            if nacids(count)~=0
                xyz_acids=[xyz_acids;xyz(count,:)];
                radius_acids=[radius_acids;sqrt(xyz(count,1)^2+xyz(count,2)^2+xyz(count,3)^2)];
            end
            imfp(count)=eventdata{count}.imfp;
            theta(count)=eventdata{count}.theta;
            phi(count)=eventdata{count}.phi;

            if plotfig~=0
                figure(fig1);
                switch act{count}
                    case 'SE'
                        plot3(eventdata{count}.xyz(1),eventdata{count}.xyz(2),eventdata{count}.xyz(3),'ob','markersize',12,'linewidth',3.0);
                    case '6eVRes'
                        plot3(eventdata{count}.xyz(1),eventdata{count}.xyz(2),eventdata{count}.xyz(3),'or','markersize',12,'linewidth',3.0);
                    case 'vibr'
                        plot3(eventdata{count}.xyz(1),eventdata{count}.xyz(2),eventdata{count}.xyz(3),'ok','markersize',12,'linewidth',3.0);
                    case 'none'
                        plot3(eventdata{count}.xyz(1),eventdata{count}.xyz(2),eventdata{count}.xyz(3),'og','markersize',12,'linewidth',3.0);
                    case '6eVRes-none'
                        plot3(eventdata{count}.xyz(1),eventdata{count}.xyz(2),eventdata{count}.xyz(3),'*r','markersize',12,'linewidth',3.0);
                end

                drawnow;
                hold on;plot3(xyz(:,1),xyz(:,2),xyz(:,3),'linewidth',3.0)
            end
        end
    end
    nacids_thrutrial(trial_count)=sum(nacids);
    tend=toc(tst_trials);
    fprintf('...Took %.2f s \n',tend);
end

tend=toc(tstart);
fprintf('Total time = %.2f s\n',tend);

fprintf('Total # of acids = %d\n',sum(nacids));



