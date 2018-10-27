%%%%----------------------Scattering Sims---------------------------%%%%
%   Note: Now it's written in a portable way that as long as one copies 
%   the ElectronInteractions completely no change is needed for the paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0 Initialization
clear;
clc;
close all;
set(0,'DefaultFigureWindowStyle','docked')

% Global functions
addpath('..\..\GlobalFunctions');

cross_sect_data_path='..\..\Traj_MonteCarlo\CrossSect_Data\';

scattdata.vibr=load([cross_sect_data_path 'VibrExcit_Data_Khakoo_2013.mat']);

scattdata.vibr.datasrc='Khakoo';

scattdata.vibr.imfp_func=@ephscatt;

optdata_path='..\..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';

scattdata.optical=load([optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData.mat']);

TruongData=load('TruongData_PMMA_PS.mat');

pathname='..\..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';

filename='DDCSdata_Ef=0p5_Elossmin=0.001eV_Erange=[5,1000].mat';

filename='DDCSdata_withICSData_Ef=10_Elossmin=0.001eV_Erange=[5,1000].mat';

filename='DDCSdata_Fuji_Ef=15.5_Elossmin=3eV_Erange=[19,200]_EQCons=Pines.mat';

filename='DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,200].mat';

scattdata.optical.inel_dcsdata=load([pathname filename]);
scattdata.E_inel_thr=min(scattdata.optical.E);

% File output base path
outputBasePath  =   strcat('..\\..\\..\\..\\JonathanCodeIO_CXRO\\',...
            'ElectronInteractions\\LEEMRes\\RealScatt_80_6_PropScatt_OffCntr\\');
        
%% Globals for tracking purposes
global  secSpawningTheta scattVector thetaLog;

%%% These initializations are necessary otherwise something might spillover
%%% even if clear command is issued (just in case)
secSpawningTheta =  [];
scattVector      =  [];
thetaLog         =  [];


% The illustration variable toggles whether explainative remarks show up
global illustration debugOutput;
illustration = 0;
debugOutput = {};
%% Debug parameters
debugCurrAcidArrayDiff = 0;
%% Model Parameters
%%% The limit where scattering ceases
SCATTERING_LOW_ENERGY_CUTOFF    =       20;
%%% The energy where the electron enters low energy regime
LOW_ENERGY_BEHAVIOUR_BOUNDARY   =       20; 
%%% Low energy random walk mean free path.
LOW_ENERGY_MEAN_FREE_PATH       =       3.67;

%% 1 Creating the PAG grid
rho_pag=0.4; % pag per nm^3
rho_polym=6; % polymer per nm^3
% rho_pag=0;
px_nm=[1 1 1]*0.25;

% univ is the volume of concern
% size of the grid in nm
univ.size_nm=[51 51 51]; % nm
% voxel size in nm
univ.px_nm=px_nm;
% size of the grid in # of voxels
univ.npx=floor(univ.size_nm./univ.px_nm);

% x,y,z coordinates of the voxels in nm
xgrid=-univ.size_nm(1)/2+px_nm(1)/2     :px_nm(1)   :univ.size_nm(1)/2-px_nm(1)/2;
ygrid=-univ.size_nm(2)/2+px_nm(2)/2     :px_nm(2)   :univ.size_nm(2)/2-px_nm(2)/2;
zgrid=-univ.size_nm(3)/2+px_nm(3)/2     :px_nm(3)   :univ.size_nm(3)/2-px_nm(3)/2;
% To generate a meshgrid whose indices correspond to axes (x,y,z)
[ygrid,xgrid,zgrid]     =       meshgrid(ygrid,xgrid,zgrid);


% converting line spaces into meshgrids
univ.grid.x=xgrid;
univ.grid.y=ygrid;
univ.grid.z=zgrid;

% pag density 1/voxel (instead of nm^3)
pagimg=rho_pag*prod(univ.px_nm).*ones(univ.npx);
% generating the pag loading by poisson randomization
pagimg=poissrnd(pagimg);

% initiating the electron image
elecimg=zeros(size(pagimg));

%% 2 Initializing the event structure [holds info about excitations triggered]
% event{1}.xyz=round(univ.npx/2);
event{1}.xyz=[0 0 0];
xyzglobal.x=event{1}.xyz(1);
xyzglobal.y=event{1}.xyz(2);
xyzglobal.z=event{1}.xyz(3);

E_in=91;

event{1}.Ein=E_in;
event{1}.Eout=0;
event{1}.Eloss=0;
% The energy of the "secondary electron" generated or the starting energy
% if a new trace is created out of this event
% The primary event is the incidence and the primary electron energy goes
% there to trigger Scarrcalc_lowE. It signifies an "energetic electron" 
% (and another event trace) is created
event{1}.Ese=91;
event{1}.act='SE';
event{1}.nacid=0;
event{1}.nSE=0;
event{1}.elecimg=elecimg;
event{1}.pag.img=pagimg;
event{1}.pag.rcnrad=3; % nm
event{1}.pag.rho=rho_pag;
% event{1}.controlparms.model_pag_devel=1;
event{1}.univ=univ;
event{1}.scatt_Elim     =   SCATTERING_LOW_ENERGY_CUTOFF;
event{1}.lowEthr        =   LOW_ENERGY_BEHAVIOUR_BOUNDARY;
event{1}.lowEimfp       =   LOW_ENERGY_MEAN_FREE_PATH;

%% 3 Scan sweep parameters
% Number of trials per energy
% ntrials=sum(absimg(:));
% abspos=find(absimg~=0);
ntrials=1000;
tstart=tic;
% Configuring Esweep
% Esweep=linspace(20,100,5);
Esweep=[80];
% Esweep=[30];
pathlen=[];
Energy=[];

%% 4 Initiating electron incidence and dose parameters
% electron incidence grid (2D)
% in pixel coordinates
elecimg_inc=zeros([univ.npx(1) univ.npx(2)]);

% The pattern of incidence
pattsize=0; % nm (5 nm by default)
% The linespaces for the incedent square
rowsel=round(univ.npx(1)/2-pattsize/2/univ.px_nm(1)):round(univ.npx(1)/2+pattsize/2/univ.px_nm(1));
colsel=round(univ.npx(2)/2-pattsize/2/univ.px_nm(2)):round(univ.npx(2)/2+pattsize/2/univ.px_nm(2));

%{
rowsel=rowsel(2:end); % remove the extra 1 addeded due to integer math
colsel=colsel(2:end);
%}

% Applying the mask onto the incidence map
elecimg_inc(rowsel,colsel)=1;
elecimg_bin=elecimg_inc;

% Dose per masking unit
% [78,157,314];
Dose=2; % /pixel;

% Convert mask into dose
elecimg_inc=Dose.*elecimg_inc;
% elecimg_inc=floor(elecimg_inc);

%%% Debug:  overrides everything above this line in this section
%%% Effect: trial sweep having only 1 electron incident at
%%% center, uncomment the 2 lines below:
%     elecimg_inc=zeros([univ.npx(1) univ.npx(2)]);
%     elecimg_inc(round(size(elecimg_inc,1)/2),round(size(elecimg_inc,2)/2))=Dose;


% Seed the random number generator
% rng(56);
rng('shuffle');
% elecimg_inc=poissrnd(elecimg_inc);
rng('shuffle');

% Display the electron inciednce dose
figure;imagesc(elecimg_inc);colorbar;title('Electron image');drawnow;

% Acquireing the coordinates of the electrons
[x_inc,y_inc]=ind2sub(size(elecimg_inc),find(elecimg_inc>=1));

%% 5 The Simulation Loops
fig1=figure;

logfile_fid=fopen('logfile.dat','w');
fprintf(logfile_fid,'Simulation started at %s\n',datestr(now));

% A loop through the energies of interest
for E_count=1:length(Esweep)
    %% 5..1 Iteration for an energy of interest
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
    acidimg_global2=zeros(univ.npx);
    ionimg_global=zeros(univ.npx);
    radius_acids_thrutrial={};
    radius_ions_thrutrial={};
    nacids_thrutrial=[];
    nacids_unsat_thrutrial=[];
    nions_thrutrial=[];
    meanacids_thrutrial=[];
    meanions_thrutrial=[];
    
    theta_init=[];
    
    % A loop through the number of trials
    for trial_count=1:ntrials
        %% 5..2 Iteration for trial (Each trial is a clean start: a new system)
        %%% Array for acid positions to the nearest voxel
        acid_xyz=[];
        %%% Array for acid activation event positions
        acid_fine_xyz=[];
        radius_acids=[];
        radius_ions=[];
        ion_xyz=[];

        tst_trials=tic;
        fprintf(logfile_fid,'\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
        fprintf('\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
        
        % The pag loading is re-initizlized here
        % And so is the polymer loading
        pagimg=rho_pag*prod(univ.px_nm).*ones(univ.npx);
        polym_img=rho_polym*prod(univ.px_nm).*ones(univ.npx);
        
        %%%% Load your own pag and polym images
%         tmpdata=load('Scratch_Folder\PAG_Polym_Images_29eVTests.mat');
%         pagimg=tmpdata.pagimg_pre;
%         polym_img=tmpdata.polymimg_orig;
        
        rng('shuffle');
        pagimg=poissrnd(pagimg);
        
        rng('shuffle');
        polym_img=poissrnd(polym_img);
        rng('shuffle');
        
        % registering the properties of the original pag loading
        pagimg_orig=pagimg;
        polymimg_orig=polym_img;
        fprintf(logfile_fid,'Initial <PAG> = %.4f/nm^3\n',mean(pagimg_orig(:))/prod(univ.px_nm));
        fprintf(logfile_fid,'Initial <Polymer> = %.4f/nm^3\n',mean(polymimg_orig(:))/prod(univ.px_nm));
        fprintf('\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
%         event{1}.pag.img=pagimg;

        % Initializing the pag info structure (pagdata)
        pagdata.pagimg=pagimg;
        pagdata.acid_act_xyz_idx    =[]; % The index of the voxel
        pagdata.acid_act_xyz        =[]; % The position of event
        % Initializing the polymer info structure (polymdata)
        polymdata.polym_img=polym_img;
        polymdata.SE_act_xyz_idx=[];
        
        % Initialize the acid image
        acidimg=zeros(size(pagimg));
        
%       xyz_acids=[];radius_acids=[];
        xyz_electron=[];
%       nacids_thrutrial=[];

        totalacids_pertrial=0;        
        nacids_unsat_total=0;
        
        % loop through all positions with incidend electrons
        for xyz_count=1:length(x_inc)
            %% 5..3 Iteration for dosed positions
            tst_coordinates=tic;
            fprintf(logfile_fid,'\n...Co-ordinate %d of %d\n',xyz_count,length(x_inc));
%             fprintf('\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
            nelectrons=elecimg_inc(x_inc(xyz_count),y_inc(xyz_count));
            
            %% Iinitial position of the electron in nm
            xval=univ.grid.x(x_inc(xyz_count),y_inc(xyz_count),1);
            yval=univ.grid.y(x_inc(xyz_count),y_inc(xyz_count),1);
%            zval=min(univ.grid.z(:)); % put the electron at the origin in z.
            zval=0; % put it in the center

            for photon_count=1:nelectrons
                %% 5..4 Iteration for a single electron
                tst_nelec=tic;
                fprintf('......Electron %d of %d\n',photon_count,nelectrons);
                fprintf(logfile_fid,'......Electron %d of %d\n',photon_count,nelectrons);
%                 fprintf('\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
                event{1}.xyz=[xval yval zval];
                
                %Suchit's initialization
                %event{1}.theta_in=-pi+pi*rand;                 
                %Jonathan's initialization
                cosTheta    =   rand(1)*2-1;
                event{1}.theta_in = acos(cosTheta);
                theta_init=[theta_init event{1}.theta_in];
                rng('shuffle');
                event{1}.phi_in=2*pi*rand; % Actually this does nothing, as in trajcalc3, first phi used is sampled from distributions
                
%               event{1}.pag.img=pagimg;
                
                xyzglobal.x=event{1}.xyz(1);
                xyzglobal.y=event{1}.xyz(2);
                xyzglobal.z=event{1}.xyz(3);
%               xyzglobal.z=0; % put the electron at the center.

                eventdata={};
                
                % Again, registering the pre scattering pag 
                % and polymer loading
                pagimg_pre=pagdata.pagimg;
                polym_img=polymdata.polym_img;
                
                % The scattering engine itself
                [eventdata,pagdata,polymdata]=Scattcalc_lowE_JHSM(event,scattdata,eventdata,xyzglobal,pagdata,polymdata,logfile_fid);

%%% These are legacy lines that come previous ways of passing pagdata %%%%%                
%               eventdata_global{trial_count,E_count}=eventdata;
%               acidimg=acidimg+pagimg-eventdata{end}.pag.img;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Rasing a flag if pag loading becomes negative or
                % increases
                if any(pagdata.pagimg(:)<0)
                    dbg=1;
                end
                acidimg=acidimg+pagimg_pre-pagdata.pagimg;
                if any(pagimg_pre(:)<pagdata.pagimg(:))
                    dbg=1;
                end
%               pagimg=eventdata{end}.pag.img;
%               pagimg(pagimg<0)=0;

                %% 5..5.1 In loop analysis of the results
                
                %%% Storing the number of acids generated in this trial
                %%% up to the current electron
                nacids_total=[nacids_total sum(acidimg(:))];

                act={};nSE=[];imfp=[];theta=[];phi=[];nacids=[];nacids_unsat_total_tmp=[];

                pathlen(trial_count,1)=0; % its 0 when electron energy is the incident energy.
                nacids_count=1;
%                 xyz_acids=[];%radius_acids=[];
                for count = 1:length(eventdata)
                    %% 4..5 Analysing each and every scattering event from
                    Ein(trial_count,count)      =eventdata{count}.Ein;
                    Eloss(trial_count,count)    =eventdata{count}.Eloss;
                    Ese(trial_count,count)      =eventdata{count}.Ese;
%                     xyz(count,:)=eventdata{count}.xyz;
                    xyz_electron                =[xyz_electron;...
                        eventdata{count}.xyz];
                    act{count}                  =eventdata{count}.act;
                    nSE(count)                  =eventdata{count}.nSE;
                    nacids(count)               =eventdata{count}.nacid;
                    nacids_unsat_total_tmp(count)...
                                                =eventdata{count}.nacid_unsat;
                    totalacids_pertrial         =totalacids_pertrial+eventdata{count}.nacid;
                    
                    pathlen(trial_count,count+1)=eventdata{count}.pathlen;
                    Energy(trial_count,count)   =eventdata{count}.Ein;
%                   Energy=[Energy eventdata{end}.Eout];
        
                    imfp(count)                 =eventdata{count}.imfp;
                    theta(count)                =eventdata{count}.theta;
                    phi(count)                  =eventdata{count}.phi;
                    
                    % Updating "global variables"
                    Energy_global       =   [Energy_global;eventdata{count}.Ein];
                    imfp_global         =   [imfp_global;eventdata{count}.imfp];
                    steplen_global      =   [steplen_global;eventdata{count}.rnew];
                    act_global{length(act_global)+1}...
                                        =   eventdata{count}.act;
                    %%% Registers the final angles of all events
                    theta_global(length(theta_global)+1)...
                                        =   eventdata{count}.theta;
                    phi_global(length(phi_global)+1)...
                                        =   eventdata{count}.phi;
                    %%% Registers the final positions of all events
                    %%% Positions are not coarse grained by pixel bin sizes
                    xyz_electron_global(size(xyz_electron_global,1)+1,1:3)...
                                        =   eventdata{count}.xyz;
                    
                end
                nacids_unsat_total=nacids_unsat_total+sum(nacids_unsat_total_tmp);
                
                Energy(trial_count,count+1)=eventdata{count}.Eout;
                nacids_total2=[nacids_total2 sum(nacids)];
                
                % tend stand for t-end of the timer
                % Recodr the time spent per elctron
                tend=toc(tst_nelec);
                fprintf(logfile_fid,'......Took %.2f s \n',tend);
            end
            
            % tend stand for t-end of the timer
            % Recodr the time spent per dosed site
            tend_coordinates=toc(tst_coordinates);
            fprintf(logfile_fid,'...Took %.2f s\n',tend_coordinates);
        end
        
        %%% unsaturated acid count
        nacids_unsat(trial_count,E_count)=nacids_unsat_total;
        
        %%% Do acid counting etc. right before next trial and energy values
        %%% are used:
        %%% The indices of the posisions where acid is generated
        acid_act_xyz_idx    =   pagdata.acid_act_xyz_idx;
        acid_act_xyz        =   pagdata.acid_act_xyz;
        %%% Update the global images
        %%% The pick-and-replace image
        acidimg_global(acid_act_xyz_idx)    =acidimg(acid_act_xyz_idx);
        %%% The culmulative image
        acidimg_global2(acid_act_xyz_idx)   =acidimg_global2(acid_act_xyz_idx)+acidimg(acid_act_xyz_idx);
        
        %%%% Logging acid positions
        %%% The acid coordinate in voxel scaling
        for acid_count=1:length(acid_act_xyz_idx)
            [actidx_x,actidx_y,actidx_z]...
                =ind2sub(size(pagdata.pagimg),acid_act_xyz_idx(acid_count));
            %%% The acid coordinate in nm scaling
            xyz_tmp=[univ.grid.x(actidx_x,actidx_y,actidx_z)...
                univ.grid.y(actidx_x,actidx_y,actidx_z)...
                univ.grid.z(actidx_x,actidx_y,actidx_z)];
            acid_xyz        =   [acid_xyz;xyz_tmp];
        end
        %%% The activation event coordinate
        acid_fine_xyz   =   [acid_fine_xyz; acid_act_xyz];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     Debug code to figure out why the fine position arrays is
        %%%     not as long as the pixelated one
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        if size(acid_xyz,1)-size(acid_fine_xyz,1)~= debugCurrAcidArrayDiff
            debugObject = {};
            debugObject.condition   =   'The difference between the two acid arrays changes';
            debugObject.level       =   'trajsim_outer';
            debugObject.coarseArrayLength   = size(acid_xyz,1);
            debugObject.fineArrayLength     = size(acid_fine_xyz,1);
            debugOutput{size(debugOutput,2)+1} = debugObject;
            debugCurrAcidArrayDiff = size(acid_xyz,1)-size(acid_fine_xyz,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        
        %%% Same thing for ions
        SE_act_xyz_idx  =polymdata.SE_act_xyz_idx;
        ionimg_global(SE_act_xyz_idx)=  ionimg_global(SE_act_xyz_idx)+1;
        for ion_count=1:length(SE_act_xyz_idx)
            [ionidx_x,ionidx_y,ionidx_z]=ind2sub(size(polymdata.polym_img),SE_act_xyz_idx(ion_count));
            xyz_tmp=[univ.grid.x(ionidx_x,ionidx_y,ionidx_z)...
                univ.grid.y(ionidx_x,ionidx_y,ionidx_z)...
                univ.grid.z(ionidx_x,ionidx_y,ionidx_z)];
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
        
        %%% this line added on 5/27/2017
        nacids_unsat_thrutrial(trial_count)=totalacids_pertrial;
        %%% end of line added on 5/27/2017
        
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
        
        %{
        %%% Default saving
        save(sprintf('LEEMRes\\Center_thetaSuchit\\Ein=%.2f_Dose=%.2fepnm2_Ef=15.5_pag-Emin=5_rcnrad=%.2f_PAG=0.4_T%d.mat',Esweep(E_count),Dose/prod(univ.px_nm(2:3)),event{1}.pag.rcnrad,trial_count));
        %}
        %%% Acid position only saving
        save(sprintf(strcat(outputBasePath,...
            'Ein=%.2f_Dose=%.2fepnm2_Ef=15.5_pag-Emin=5_rcnrad=%.2f_PAG=0.4_T%d.mat'),...
            Esweep(E_count),Dose/prod(univ.px_nm(2:3)),event{1}.pag.rcnrad,trial_count),...
        'acid_xyz','acid_fine_xyz','xyz_electron_global','ion_xyz');

%         tend=toc(tst_trials);
%         fprintf('...Took %.2f s \n',tend);
    end
    
    meanAcids_thruE(E_count)=mean(nacids_thrutrial);
    stdAcids_thruE(E_count)=std(nacids_thrutrial);
    
%     save(sprintf('QE_Stats_Sims_7\\QEStats_Dose=%.2f_RcnRad=%.2f_PAG=%.2f_E=%.2feV_pagEa=5_Scatt-Elim=20_PAG-EaMin=5_%dTrials_RcnThenMove.mat',Dose,event{1}.pag.rcnrad,rho_pag,Esweep(E_count),ntrials));

    figure(fig1);imagesc(squeeze(sum(acidimg_global,2)));title('sum(acidimg global,2)'); drawnow;
    
    tend=toc(tstart);
    fprintf(logfile_fid,'Total time = %.4f s = %.4f min = %.4f hours\n',tend,tend/60,tend/3600);
    fprintf('Total time = %.4f s = %.4f min = %.4f hours\n',tend,tend/60,tend/3600);

    fprintf(logfile_fid,'Total # of acids = %d\n',sum(nacids));
    fprintf('Total # of acids = %d\n',sum(nacids));

    fprintf(logfile_fid,'%d; <Acids> = %.4f; sig_Acids = %.4f\n',ntrials,mean(nacids_thrutrial),std(nacids_thrutrial));
    fprintf(1,'%d trials; <Acids> = %.4f; sig_Acids = %.4f\n',ntrials,mean(nacids_thrutrial),std(nacids_thrutrial));
end

figure;
plot(Esweep,meanAcids_thruE,'-o','linewidth',3.0,'markersize',12);
% errorbar(Esweep,meanAcids_thruE,stdAcids_thruE./sqrt(1200),'-o','linewidth',3.0,'markersize',12);
xlabel('E (eV)');
ylabel('Mean number of acids');
set(gca,'fontsize',30,'linewidth',3.0);

save(sprintf(strcat(outputBasePath,...
            'WorkSpace')));