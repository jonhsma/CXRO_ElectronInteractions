%%%%----------------------Scattering Sims---------------------------%%%%
%   Note: Now it's written in a portable way that as long as one copies 
%   the ElectronInteractions completely no change is needed for the paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%---------------------Notes of _NoCoarseGrain_ Variant------------------
%%%% This version involves a big change in data structure. No coarse
%%%% graining is involved in the simulation itself. Few new functions are
%%%% written to make the code more modulari so using it with previous
%%%% versions of functions might not be possible

%% 0.0 Path and other basic initialization
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

%% 0.0.1 File output base path
outputBasePath  =   strcat('..\\..\\..\\..\\JonathanCodeIO_CXRO\\',...
            'ElectronInteractions\\LEEMRes\\NoCoarseGrainTest_greasePan\\');
        
%% 0.1 Globals for tracking
global  secSpawningTheta scattVector thetaLog;

%%% These initializations are necessary otherwise something might spillover
%%% even if clear command is issued (just in case)
secSpawningTheta =  [];
scattVector      =  [];
thetaLog         =  [];

%% 0.2 Illustrations and Echo
% The illustration variable toggles whether explainative remarks show up
global illustration echoConfig;
illustration = 0;

%%% Configuration for echos
echoConfig.acid.perTrial    =   1;
echoConfig.acid.perElectron =   1;
echoConfig.acid.perTraj     =   1;

%%% Default graph numbers
FIGURE_TRAJ_PER_TRIAL       =   7201;
FIGURE_TRAJ_PER_ENERGY      =   7202;
%% 0.3 Debug parameters
global debugOutput ;
debugOutput = {};
debugCurrAcidArrayDiff = 0;
%% 1.0 Model Parameters
%%% The limit where scattering ceases
SCATTERING_LOW_ENERGY_CUTOFF    =       20;
%%% The energy where the electron enters low energy regime
LOW_ENERGY_BEHAVIOUR_BOUNDARY   =       20; 
%%% Low energy random walk mean free path.
LOW_ENERGY_MEAN_FREE_PATH       =       3.67;
%%% The reaction radius of PAGS
ACID_REACTION_RADIUS            =       3;
%% 1.1 System specification
%%% These are constants so this is the only time where they are on the LHS
%%% The x,y,z sizes of the system in nm
SYSTEM_SIZE         =   [51; 51; 51];
%%% Zero the center of the grid and set the coundaries
SYSTEM_LIMITS(:,1)  =   -SYSTEM_SIZE/2;
SYSTEM_LIMITS(:,2)  =   SYSTEM_SIZE/2;
%%% Density of photoacid generator (per nanometer cube)
RHO_PAG         =   0.4;
%%% Density of polymer  (per nanometer cube)
RHO_POLYMER     =   6;
%% 1.2 Initializing the event structure [holds info about excitations triggered]
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
%event{1}.elecimg=elecimg;
%event{1}.pag.img=pagimg;
event{1}.pag.rcnrad     =   ACID_REACTION_RADIUS; % nm
%event{1}.pag.rho=rho_pag;
%event{1}.controlparms.model_pag_devel=1;
%event{1}.univ=univ;
event{1}.scatt_Elim     =   SCATTERING_LOW_ENERGY_CUTOFF;
event{1}.lowEthr        =   LOW_ENERGY_BEHAVIOUR_BOUNDARY;
event{1}.lowEimfp       =   LOW_ENERGY_MEAN_FREE_PATH;
%% 2 Scan sweep parameters
% Number of trials per energy
nTrials=10;
eSweep=[80];

tStart=tic;

pathlen=[];
Energy=[];
%% 3.1 Initiating electron incidence and dose parameters
%%% No-matter-what-you're-using parameters
nElectrons = 10; % number of electron per trial

%% 3.1.1 Stochastic volumetric dosing parameters
dosingLimits =... The space within which dosing occurs
    [-0.5,-0.5,-0.5;...
    0.5,0.5,0.5]';
% Is the total number of electrons definite (1) or a Poisson number (0)
absoluteDosing = 0;

%% 3.1.2 Dosing trajectory generator
% The dosing trajectory (1:k) is the x position of the k-th electron
dosingPath = zeros([3 100]);
    
% The dose number that will show up in the file names
manualDose = 0;
% If array is not empty the path will be used automatically unless mandate
% by mandatoryRandom (which allows you to keep the path)
mandatorySequenceRandom = 0;
%% 3.2 Generatlized electron dosing sequence generator
% The dosingSequence Handle will give an array of size 3-by-n nomatter
% which method is selected for each trial. Function handle allows random distribution to
% be randomized for every trial
if (size(dosingPath,1)==3) && ~mandatorySequenceRandom
    dosingSequenceHandle = @(n)dosingPath(:,1:min(n,size(dosingPath,2)));
    dose = manualDose;
elseif absoluteDosing == 1
    dosingSequenceHandle = @(n)randPosGen(dosingLimits,-n);
    dose = -nElectrons/prod(dosingLimits(:,2)-dosingLimits(:,1));
else
    dosingSequenceHandle = @(n)randPosGen(dosingLimits,n);
    dose = nElectrons/prod(dosingLimits(:,2)-dosingLimits(:,1));
end
%% 4 The Simulation Loops
%fig1=figure;
logfile_fid=fopen('logfile.dat','w');
fprintf(logfile_fid,'Simulation started at %s\n',datestr(now));

%%% Results storage array
scanArchive = cell([length(eSweep) nTrials]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Energy Resolved Registers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanAcids_thruE     =   zeros([1 length(eSweep)]);
stdAcids_thruE      =   zeros([1 length(eSweep)]);

% A loop through the energies of interest
for E_count=1:length(eSweep)
    %% 4..1 Iteration for an energy of interest
    event{1}.Ein=eSweep(E_count);
    event{1}.Ese=eSweep(E_count);
        
    eventdata_global={};
    Eloss=[];

    nAcidsTotal    =  [];
    nacids_total2=[];
    Energy_global=[];
    imfp_global=[];
    steplen_global=[];
    act_global={};
    theta_global=[];
    phi_global=[];
    xyz_electron_global=[];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Trial Resolved Registers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Acid registers
    %%% The number of acid generated in each trial
    nAcids_thruTrial    =   zeros([1 nTrials]);
    %%% The average number of acids generated by one electron. created to
    %%% accomodate stocastic dosing
    meanAcids_thruTrial    =    zeros([1 nTrials]);
    %%% The distance of each acid from origin
    radius_acids_thrutrial =    cell([1 nTrials]);
    
    
    %%% Ion registers
    nIons_thruTrial         =   zeros([1 nTrials]);
    radius_ions_thrutrial   =   cell([1 nTrials]);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Trial-acculmulative arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta_init              =   [];
    radius_acids_global     =   [];
    radius_ions_global      =   [];
    posAcid_TrAccu          =   [];
    posAcidAct_TrAccu       =   [];
    
    %% 4..2 Loop through the number of trials
    for trial_count=1:nTrials
        %% 4..2.1 Iteration for trial (Each trial is a clean start: a new system)
        %%% Array for acid positions of the acid
        acid_xyz=[];
        %%% Array for acid activation event positions
        acid_e_xyz=[];
        
        %%% Still need to figure out
        radius_acids=[];
        radius_ions=[];
        ion_xyz=[];

        tst_trials=tic;
        fprintf(logfile_fid,'\nEnergy %d of %d; Trial %d of %d\n',E_count,length(eSweep),trial_count,nTrials);
        fprintf('\nEnergy %d of %d; Trial %d of %d\n',E_count,length(eSweep),trial_count,nTrials);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Legacy pag and polymer distribution generation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
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
        %}      
        %% 4..2.2 Pag and polymer positions generation
        
        posPAG      = randPosGen(SYSTEM_LIMITS,RHO_PAG*prod(SYSTEM_SIZE));
        posPolymer  = randPosGen(SYSTEM_LIMITS,RHO_POLYMER*prod(SYSTEM_SIZE));
        
        %%% registering the properties of the original pag loading
        %%% Commented out as the notion of a grid is gone. Need a new way
        %%% to give the same information
        %{
        pagimg_orig=pagimg;
        polymimg_orig=polym_img;
        
        fprintf(logfile_fid,'Initial <PAG> = %.4f/nm^3\n',mean(pagimg_orig(:))/prod(univ.px_nm));
        fprintf(logfile_fid,'Initial <Polymer> = %.4f/nm^3\n',mean(polymimg_orig(:))/prod(univ.px_nm));
        fprintf('\nEnergy %d of %d; Trial %d of %d\n',E_count,length(Esweep),trial_count,ntrials);
        event{1}.pag.img=pagimg;
        
        %pagdata.pagimg              =   pagimg;
        %polymdata.polym_img         =   polym_img;
        
        % Initialize the acid image
        %acidimg=zeros(size(pagimg));
        
        %xyz_acids=[];radius_acids=[];
        %nacids_thrutrial=[];
        %}

        % Initializing the pag info structure (pagdata)
        
        pagdata.posPAG              =   posPAG;
        pagdata.posPAG_removed      =   [];
        pagdata.acid_act_xyz_idx    =   []; % The index of the voxel
        pagdata.acid_act_e_xyz      =   []; % The position of event
        
        % Initializing the polymer info structure (polymdata)
        
        polymdata.posPolymer        =   posPolymer;
        polymdata.SE_act_xyz        =   [];
        
        xyz_electron=[];

        totalacids_pertrial=0;        
        nacids_unsat_total=0;        
        %% 4..3 Loop through all incidence of electrons
        
        %%% The sequence of positions to be dosed
        dosingSequence     =   dosingSequenceHandle(nElectrons);
                
        %%% Initializing the archive element (delayed till this point to
        %%% accomodate stochastic dosing (different number of electron per
        %%% each trial
        scanArchive{trial_count,E_count}    =   cell([1 size(dosingSequence,2)]);
        
        for incidence   =   1:size(dosingSequence,2)
            %% 4..3.1 Iteration for dosed positions
            tst_incidence=tic;
            fprintf(logfile_fid,'\n...Co-ordinate %d of %d\n',...
                incidence,size(dosingSequence,2));

            xval = dosingSequence(1,incidence);
            yval = dosingSequence(2,incidence);
            zval = dosingSequence(3,incidence);

            fprintf('......Electron %d of %d\n',incidence,nElectrons);
            fprintf(logfile_fid,'......Electron %d of %d\n',incidence,nElectrons);
            event{1}.xyz=[xval yval zval];

            %Suchit's initialization
            %event{1}.theta_in=-pi+pi*rand;                 
            %Jonathan's initialization
            cosTheta    =   rand(1)*2-1;
            event{1}.theta_in = acos(cosTheta);
            theta_init  =   [theta_init event{1}.theta_in];
            rng('shuffle');
            event{1}.phi_in=2*pi*rand; % Actually this does nothing, as in trajcalc3, first phi used is sampled from distributions

            %%% Initialize the array that stores all events from this
            %%% electron
            eventdata={};

            % Again, registering the pre scattering pag 
            % and polymer loading
            %pagimg_pre      =   pagdata.pagimg;
            posPAG_init     =   pagdata.posPAG;
            posPAG_rmv_init =   pagdata.posPAG_removed;
            %polym_img       =   polymdata.polym_img;
            posPolymer_init =   polymdata.posPolymer;

            % The scattering engine itself
            [eventdata,pagdata,polymdata]=TrajectoryWrapper(event,scattdata,eventdata,xyzglobal,pagdata,polymdata,logfile_fid);

            % Display the location of acids created by this electron.
            if echoConfig.acid.perElectron
                fprintf('Positions of acids activted by this electron\n');
                disp(pagdata.posPAG_removed(:,...
                    size(posPAG_rmv_init,2)+1:...
                    size(pagdata.posPAG_removed,2))');
            end
            %% 4..3.1 Documenting the results
            
            scanArchive{E_count,trial_count}{incidence} = eventdata;

            % Record the time spent on each incidence
            tend_coordinates=toc(tst_incidence);
            fprintf(logfile_fid,'...Took %.2f s\n',tend_coordinates);
        end
        %% 4..2.3 Per Trial Logging and counting
        
        %%% Do acid counting etc. right before next trial and energy values
        %%% are used:
        %%% The indices of the posisions where acid is generated
        acid_act_xyz_idx    =   pagdata.acid_act_xyz_idx;
        %Tranpose for compatibility with comparison code
        acid_xyz            =   pagdata.posPAG_removed';
        acid_act_e_xyz      =   pagdata.acid_act_e_xyz;
        %%% For consistency with legacy naming
        acid_fine_xyz       =   acid_act_e_xyz;
        
        posAcid_TrAccu      =   [posAcid_TrAccu; acid_xyz];
        posAcidAct_TrAccu   =   [posAcidAct_TrAccu; acid_fine_xyz];
        
        %%% Acid counting for the trial
        nAcids_thruTrial(trial_count)       =   size(acid_xyz,1);
        meanAcids_thruTrial(trial_count)    =   size(acid_xyz,1)/size(dosingSequence,2);
        
        %%% Update the global images
        %{
        %%% The pick-and-replace image
        acidimg_global(acid_act_xyz_idx)    =acidimg(acid_act_xyz_idx);
        %%% The culmulative image
        acidimg_global2(acid_act_xyz_idx)   =acidimg_global2(acid_act_xyz_idx)+acidimg(acid_act_xyz_idx);
        %}
        
        %%% Acid echo
        if echoConfig.acid.perTrial
            fprintf('Positions of acids activted in this trial (acid_xyz)\n');
            disp(acid_xyz)
            fprintf('Positions of acids activted in this trial (NaN counting)\n');
            disp(posPAG(:,isnan(pagdata.posPAG(1,:)))')
            fprintf('Positions of acids activted in this trial (array_indexing)\n');
            disp(posPAG(:,acid_act_xyz_idx)')
        end
        
        
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
        SE_act_xyz  =   polymdata.SE_act_xyz;
        ion_xyz=[ion_xyz;SE_act_xyz'];
        
        nIons_thruTrial(trial_count)        =   size(SE_act_xyz,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Legacy radius statistics
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Warning: In these ad-hoc statistics, the radii are calculated
        %%% relative to xval, yval and zval, which are the coordinates of
        %%% the last incidence. So unless all the incidences happen at the
        %%% same spot, these numbers don't mean much
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nIons_thruTrial(trial_count)>0
            radius_ions=sqrt((ion_xyz(:,1)-xval).^2+...
                (ion_xyz(:,2)-yval).^2+...
                (ion_xyz(:,3)-zval).^2);
            radius_ions_thrutrial{trial_count}=radius_ions;
            radius_ions_global=[radius_ions_global;radius_ions];
        else
            radius_ions_thrutrial{trial_count}=NaN;
        end       
    
        
        if nAcids_thruTrial(trial_count)>0
            radius_acids=sqrt((acid_xyz(:,1)-xval).^2+...
                (acid_xyz(:,2)-yval).^2+...
                (acid_xyz(:,3)-zval).^2);
            radius_acids_thrutrial{trial_count}=radius_acids;
            radius_acids_global=[radius_acids_global,radius_acids'];
        else
            radius_acids_thrutrial{trial_count}=NaN;
        end
        dbg=1;
        
        %%% Selective Saving
        save(sprintf(strcat(outputBasePath,...
            'Ein=%.2f_Dose=%.2fepnm2_Ef=15.5_pag-Emin=5_rcnrad=%.2f_PAG=0.4_T%d.mat'),...
            eSweep(E_count),dose,event{1}.pag.rcnrad,trial_count),...
            'acid_xyz','acid_fine_xyz',...
            'xyz_electron_global',...
            'ion_xyz');            
        %% 4..2.4 Per Trial Trajectory Echo
        figure(FIGURE_TRAJ_PER_TRIAL);
        grid on
        hold off
        lengenEntries = cell(size(dosingSequence,2)+1);
        for inIterator = 1:size(dosingSequence,2)
            qX = [];
            qY = [];
            qZ = [];
            qU = [];
            qV = [];
            qW = [];
            nEvents     =   length(...
                scanArchive{E_count,trial_count}{inIterator});
            for eventIterator = 1:nEvents
                this = scanArchive{E_count,trial_count}{inIterator}{eventIterator};
                qX = [qX this.xyz_init(1)];
                qY = [qY this.xyz_init(2)];
                qZ = [qZ this.xyz_init(3)];
                qU = [qU this.xyz(1)];
                qV = [qV this.xyz(2)];
                qW = [qW this.xyz(3)];
            end
            quiver3(qX,qY,qZ,qU-qX,qV-qY,qW-qZ,0,'-','LineWidth',2);
            hold on
            legendEntries{inIterator} = strcat('Electron ',num2str(inIterator));
        end
        if length(acid_xyz) >=1
            quiver3(acid_act_e_xyz(:,1),acid_act_e_xyz(:,2),acid_act_e_xyz(:,3),...
                acid_xyz(:,1) - acid_act_e_xyz(:,1),...
                acid_xyz(:,2) - acid_act_e_xyz(:,2),...
                acid_xyz(:,3) - acid_act_e_xyz(:,3),...
                0,'--');
            scatter3(acid_xyz(:,1),acid_xyz(:,2),acid_xyz(:,3),'o')            
            legendEntries{size(dosingSequence,2)+1} = 'Acid Activation Moment';
        end
        legend(legendEntries);
        title('Trajectories and activation event in the last iteration')
        daspect([1 1 1])
        drawnow;
    end
    
    meanAcids_thruE(E_count)    =   mean(nAcids_thruTrial);
    stdAcids_thruE(E_count)     =   std(nAcids_thruTrial);
    
    tend=toc(tStart);
    fprintf(logfile_fid,'Total time = %.4f s = %.4f min = %.4f hours\n',tend,tend/60,tend/3600);
    fprintf('Total time = %.4f s = %.4f min = %.4f hours\n',tend,tend/60,tend/3600);

    fprintf(logfile_fid,'%d; <Acids> = %.4f; sig_Acids = %.4f\n',nTrials,mean(nAcids_thruTrial),std(nAcids_thruTrial));
    fprintf(1,'%d trials; <Acids> = %.4f; sig_Acids = %.4f\n',nTrials,mean(nAcids_thruTrial),std(nAcids_thruTrial));    
    %% 4..1.3 Per Energy Trajectory Echo
    figure(FIGURE_TRAJ_PER_ENERGY);
    grid on
    hold off
    qX = [];
    qY = [];
    qZ = [];
    qU = [];
    qV = [];
    qW = [];
    for trialIter   = 1:nTrials
        for inIterator = 1:size(scanArchive{E_count,trialIter},2)
            nEvents     =   length(...
                scanArchive{E_count,trialIter}{inIterator});
            for eventIterator = 1:nEvents
                this = scanArchive{E_count,trialIter}{inIterator}{eventIterator};
                qX = [qX this.xyz_init(1)];
                qY = [qY this.xyz_init(2)];
                qZ = [qZ this.xyz_init(3)];
                qU = [qU this.xyz(1)];
                qV = [qV this.xyz(2)];
                qW = [qW this.xyz(3)];
            end
        end
    end
    quiver3(qX,qY,qZ,qU-qX,qV-qY,qW-qZ,0,'.-','LineWidth',1);
    hold on            
    if length(posAcid_TrAccu) >=1
        quiver3(posAcidAct_TrAccu(:,1),posAcidAct_TrAccu(:,2),posAcidAct_TrAccu(:,3),...
            posAcid_TrAccu(:,1) - posAcidAct_TrAccu(:,1),...
            posAcid_TrAccu(:,2) - posAcidAct_TrAccu(:,2),...
            posAcid_TrAccu(:,3) - posAcidAct_TrAccu(:,3),...
            0,'--');
        scatter3(posAcid_TrAccu(:,1),posAcid_TrAccu(:,2),posAcid_TrAccu(:,3),...
            'o')
    end
    title({strcat('Acid and Trajectory visuallization at h\nu = ',num2str(eSweep(E_count)),' eV');...
        strcat(num2str(nTrials),' trials with',' ',...
        num2str(nElectrons),' electrons per trial on average')});
    legend({'Trajectories','Activations','Acids'});
    daspect([1 1 1]);
    drawnow;
    
end
%% Final workspace export
%figure;
%plot(eSweep,meanAcids_thruE,'-o','linewidth',3.0,'markersize',12);
% errorbar(Esweep,meanAcids_thruE,stdAcids_thruE./sqrt(1200),'-o','linewidth',3.0,'markersize',12);
%xlabel('E (eV)');
%ylabel('Mean number of acids');
%set(gca,'fontsize',30,'linewidth',3.0);

save(sprintf(strcat(outputBasePath,...
            'WorkSpace')));