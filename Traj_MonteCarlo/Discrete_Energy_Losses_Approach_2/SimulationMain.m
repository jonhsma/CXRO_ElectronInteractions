%%%%---------------Scattering Trajectory Simulation-------------------%%%%
% This is a Monte Carlo electron simulation. Electrons are generated and
% the evolution (propagation, scattering and termination) of individual
% electrons are simulated in a step-by-step manner. The program is
% organized is the following way. This script configures the system,
% initiate simulations and handles the data logistics. Each initiation is
% sent to the trajectory wrapper. The wrapper does two things. It sends the
% initiation to the trajectory follower, which simulates the life of the
% electron (which is pretty short). When the simulation is over, the
% wrapper reads in the trajectory, goes through each step and sees if a
% secondary is generated. If so, it does it first job again and initiate
% another trajectory. The trajectory follower is where the
% physics/chemistry is simulated. It simulates the life of the elctron and
% returns it's trajectory and relavant information.
%=========================================================================
% Major changes to the flow and structure in chronological order
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Code was inherited from Suchit -SB
%-------------------------------------------------------------------------
% Absolute file depedence is abolished. Should run if you move the
% ElectronInteractions folder as a whole - JHM
%-------------------------------------------------------------------------
% Voxelization of acid/polymer loading is abandoned for sub-nm resolution.
% No side effect on performance observed - JHM
%-------------------------------------------------------------------------
% Flexible dosing and angular distribution. Both of them are configurable
% as they are now function handles. -JHM
%------------------------------------------------------------------------- 
% Parallelization of cycles. When debugging, turn off parallelization. New
% graphical tools developed to observe trajectories after the fact. - JHM
%-------------------------------------------------------------------------
% Use of optional trajectory follower allowed. JHM
%=========================================================================
% Architectural work in progress
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Cleaning up and reorganizing the repository
%-------------------------------------------------------------------------
% Make scattering engines more plug-and-play in the sense that changes
% needed to be made to the follower is minimal when toggling a scattering
% channel
%=========================================================================

%% 0.0 Path and other basic initialization
clear;
clear global
clc;
close all;
set(0,'DefaultFigureWindowStyle','docked')


% Global functions
addpath('..\..\GlobalFunctions');

% Components
addpath('..\CommonComponents');

cross_sect_data_path='..\..\Traj_MonteCarlo\CrossSect_Data\';

scattdata.vibr=load([cross_sect_data_path 'Vibr_Khakoo_2013_x50.mat']);
%%% Determining the vibrational chanel of scattering
scattdata.vibr.datasrc='Khakoo';
%%% The Frolich IMFP function
scattdata.vibr.imfp_func=@ephscatt;
%%% Optical data path
optdata_path='..\..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';

scattdata.optical=load([optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData.mat']);

TruongData=load('TruongData_PMMA_PS.mat');

pathname='..\..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';
%{
filename='DDCSdata_Ef=0p5_Elossmin=0.001eV_Erange=[5,1000].mat';

filename='DDCSdata_withICSData_Ef=10_Elossmin=0.001eV_Erange=[5,1000].mat';

filename='DDCSdata_Fuji_Ef=15.5_Elossmin=3eV_Erange=[19,200]_EQCons=Pines.mat';
%}
filename='DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,200].mat';

scattdata.optical.inel_dcsdata  = load([pathname filename]);
scattdata.E_inel_thr            = min(scattdata.optical.E);

%% 0.0.1 --> File output base path
outputParent    =   strcat('..\\..\\..\\..\\JonathanCodeIO_CXRO\\',...
            'ElectronInteractions\\LEEMRes\\');
outputFolder    =   '20190111_PES_FixedIMFP';
outputBasePath  =   strcat(outputParent,outputFolder,'\\');

%% 0.0.2 --> Scattering Engines' Paths
%%% The subfolders that stores the scattering codes
scattEnginePaths={...
    'OptDataScatt\';...
    'VibrationScatt\';...
    'StoneWall\'};

%% 0.0.3 Folder management
for iPath = 1:size(scattEnginePaths)
    addpath(scattEnginePaths{iPath})
end

outputBasePathConfirmed = 0;

while ~outputBasePathConfirmed
    if ~exist(outputBasePath,'file')
           mkdir(outputBasePath);
           outputBasePathConfirmed  =   1;
    else
        prompt = strcat('Folder already exists.',...
            'Please confirm the folder name.',...
            'WARNING: FILES IN THE EXISTING FOLDER WILL BE OVERWRITTEN. ',...
            'Press cancel to abort execution');
        promptTitle = 'Possible output data folder conflict';
        prmptDims = [1 120];
        promptDefInput = {outputFolder};
        promptAnswer = inputdlg(prompt,promptTitle,prmptDims,promptDefInput);

        if max(size(promptAnswer))<1
            return;
        end
        if strcmp(outputFolder,promptAnswer{1})
            outputBasePath  =   strcat(outputParent,outputFolder,'\\');
            mkdir(outputBasePath);
            outputBasePathConfirmed = 1;
        else 
            outputFolder = promptAnswer{1};
            outputBasePath  =   strcat(outputParent,outputFolder,'\\');
        end
    end
end
diary(strcat(outputBasePath,'ConsoleLog'));
diary on;

%% 0.1 Globals for tracking
global  secSpawningTheta scattVector thetaLog;

%%% These initializations are necessary otherwise something might spillover
%%% even if clear command is issued (just in case)
secSpawningTheta =  [];
scattVector      =  [];
thetaLog         =  [];

%% 0.2 --> Illustrations and Echo
% The illustration variable toggles whether explainative remarks show up
% In Parallelization, perTrial echos are not possible
global illustration echoConfig;
illustration = 0;

%%% Configuration for echos
echoConfig.acid.perTrial    =   0;
echoConfig.acid.perElectron =   0;
echoConfig.acid.perTraj     =   0;

echoConfig.traj3.perTrial   =   0;
echoConfig.traj3.perEnergy  =   0;

echoConfig.acidDist.active  =   1;
echoConfig.acidDist.range   =   20; % 10 means from -10 to 10 nm,
echoConfig.acidDist.res     =   1;

%%% Default graph numbers
FIGURE_TRAJ_PER_TRIAL       =   7201;
FIGURE_TRAJ_PER_ENERGY      =   72000;
FIGURE_ACID_DIST_PERENERGY  =   72100;
FIGURE_ACID_VS_ENERGY       =   1001;

%% 0.3 Debug Grease Pan
global debugOutput ;
debugOutput = {};
debugCurrAcidArrayDiff = 0;

%% 1.0 --> Model Parameters
%--------------------------------------------------------------------------
% TRAJECTORY FOLLOWER OPTIONS
%--------------------------------------------------------------------------
%%% The trajectory follower options. Anything that's not a proper function
%%% handle will result in using the default trajectory follower
scattdata.followerHandle        =       @TrajectoryFollowerFixedIMFP;%'default';@TrajectoryFollowerRandomWalk@TrajectoryFollowerLowEnergyAcid;;
%--------------------------------------------------------------------------

% GENERAL SCATTERING OPTIONS
%%% The limit where scattering ceases
SCATTERING_LOW_ENERGY_CUTOFF    =       2;
%%% The energy where the electron enters low energy regime
%%% (energy below which the optical MFP is reaplced by the constant LOW_ENERGY_MEAN_FREE_PATH)
% Now defunct. DOESN'T DO ANYTHING
LOW_ENERGY_BEHAVIOUR_BOUNDARY   =       20; 
%%% Low energy mean free path of plasmon scattering.
LOW_ENERGY_MEAN_FREE_PATH       =       3.67;%1 for random walk test;
%%% The reaction radius of PAGS
ACID_REACTION_RADIUS            =       3;

% VIBRATIONAL OPTIONS
%%% Molecular density for vicrational calculations
MOLECULAR_NUMBER_DENSITY        =       1.2/120*6.02*1e23*10; % molecules/cm3

% STONEWALL OPTIONS
%%%% The energy below which the electron would activate an acid and die
scattdata.stoneWall.CUTOFF      =   5;
%%%% The imfp of the stone wall. Should be small if active
scattdata.stoneWall.IMFP        =   0.001;  
%%%% The reaction radius of the stone wall
scattdata.stoneWall.ACID_REACTION_RADIUS   =   3;  

%% 1.1 --> System specification
%%% These are constants so this is the only time where they are on the LHS
%%% The x,y,z sizes of the system in nm
SYSTEM_SIZE         =   [51; 51; 51];
%%% Zero the center of the grid and set the coundaries
SYSTEM_LIMITS(:,1)  =   -[SYSTEM_SIZE(1:2)/2; SYSTEM_SIZE(3)];
SYSTEM_LIMITS(:,2)  =   [SYSTEM_SIZE(1:2)/2; 0];
%SYSTEM_LIMITS = [-SYSTEM_SIZE/2, SYSTEM_SIZE/2];
%%% Density of photoacid generator (per nanometer cube)
RHO_PAG         =   0.4;
%%% Density of polymer  (per nanometer cube)
RHO_POLYMER     =   6;
%%% Pass it into the trajectory follower
scattdata.SYSTEM_LIMITS     =   SYSTEM_LIMITS;

%% 1.2 Initializing the event structure [holds info about excitations triggered]
% The energy of the "secondary electron" generated or the starting energy
% if a new trace is created out of this event
% The primary event is the incidence and the primary electron energy goes
% there to trigger Scarrcalc_lowE. It signifies an "energetic electron" 
% (and another event trace) is created

eventPrototype.xyz=[0 0 0];
xyzglobal.x=eventPrototype.xyz(1);
xyzglobal.y=eventPrototype.xyz(2);
xyzglobal.z=eventPrototype.xyz(3);

E_in=91;

eventPrototype.Ein=E_in;
eventPrototype.Eout=0;
eventPrototype.Eloss=0;
eventPrototype.Ese=91;
eventPrototype.act='SE';
eventPrototype.nacid=0;
eventPrototype.nSE=0;
eventPrototype.pag.rcnrad     =   ACID_REACTION_RADIUS; % nm
eventPrototype.scatt_Elim     =   SCATTERING_LOW_ENERGY_CUTOFF;
eventPrototype.lowEthr        =   LOW_ENERGY_BEHAVIOUR_BOUNDARY;
eventPrototype.lowEimfp       =   LOW_ENERGY_MEAN_FREE_PATH;

%% 2   --> Scan sweep parameters
% Number of trials per energy
nTrials     =   176;
eSweep      =   [85];% 80 65 50 40 30];
tStart      =   tic;

%% 3.1 --> Initial electron incidence and dose parameters
%%% No-matter-what-you're-using parameters
nElectrons = 100; % number of electron per trial

%% 3.1.1 --> Stochastic volumetric dosing parameters
dosingLimits =... The space within which dosing occurs
    [-2.5,2.5;...
    -2.5,2.5;...
    -10,0.0];
% Is the total number of electrons definite (1) or a Poisson number (0)
absoluteDosing = 1;

%% 3.1.2 --> Dosing trajectory generator
% The dosing trajectory (1:k) is the x position of the k-th electron
dosingPath = zeros([3 100]);
    
% The dose number that will show up in the file names
manualDose = 1;
% If array is not empty the path will be used automatically unless mandate
% by mandatoryRandom (which allows you to keep the path)
mandatorySequenceRandom = 1;

%% 3.1.3 Generatlized electron dosing sequence generator
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

%% 3.2 Electron angular distribution handle 
% One can specify the per unit solid angle distribution (as a function of
% theta) of the incident electrons

%%% Resolution of theta axis. The PDF is generated once so the resolution
%%% of the theta axis shouldn't matter
incTGThetaN    = 1000;   
incTGThetaAxis = 0:pi/incTGThetaN:pi;

%% 3.2.1 --> The probability distribution. Need not worry about normalization
%%% DO NOT sine weight as sine weighing is included in the CDF generation
incTGPDF       = zeros([1 size(incTGThetaAxis,2)]);%sin(incTGThetaAxis).^38;
incTGPDF(:)    = 1; %incTGPDF(size(incTGPDF,2)) = 1;
%%% This distribution correspond to a SEM style downward injection
%% 3.2.2 Obtain the culmulative ditribution
incTGCDF       = zeros([1 incTGThetaN+1]);
for i = 2:incTGThetaN+1
    incTGCDF(i)       = trapz(incTGThetaAxis(1:i),...
        sin(incTGThetaAxis(1:i)).*incTGPDF(1:i));
end
%%% Ensur the culmulative distribution goes to 1 at pi
incTGCDF    =   incTGCDF/incTGCDF(end);
%%% Declear the handle
thetaGenHandle = @(x)randgen(incTGThetaAxis,incTGCDF,x); 

%% 4 The Simulation Loops

logfile_fid=fopen('logfile.dat','w');
fprintf(logfile_fid,'Simulation started at %s\n',datestr(now));

%%% Results storage array
scanArchive = cell([1 length(eSweep)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Energy Resolved Registers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanAcids_thruE     =   zeros([1 length(eSweep)]);
stdAcids_thruE      =   zeros([1 length(eSweep)]);

% Record the workspace prio to the the simulation
save(sprintf(strcat(outputBasePath,...
            'preExecutionWS')));

%% 4.1 Loop through the energies of interest
for E_count=1:length(eSweep)
    %% 4..1 Iteration for an energy of interest
    eventPrototype.Ein=eSweep(E_count);
    eventPrototype.Ese=eSweep(E_count);
        
    eventdata_global={};
    Eloss=[];

    nAcidsTotal    =  [];
    nacids_total2=[];
    Energy_global=[];
    imfp_global=[];
    stepLen_global=[];
    act_global={};
    theta_global=[];
    phi_global=[];
    xyz_electron_global=[];
    
    %%% Timing
    tStart_E = tic;
    
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
    ionXyz_TrAccu           =   [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Parallelism adaptation
    energyScanArchive   =   cell([1 nTrials]);    
    
    %% 4..2 Loop through the number of trials
    fprintf('\n Simulations with %4.2f eV incident electrons started\n',...
        eSweep(E_count));
    parfor trial_count=1:nTrials
    %for trial_count=1:nTrials
        %% 4..2.1 Iteration for trial (Each trial is a clean start: a new system)
        
        %%% Logging metadata
        energyScanArchive{trial_count}.energy         = eSweep(E_count);
                
        %%% Array for acid positions of the acid
        acid_xyz    =   [];
        %%% Array for acid activation event positions
        acid_e_xyz  =   [];
        
        %%% Still need to figure out
        radius_acids=   [];
        radius_ions=[];
        ioyn_xz=[];

        tst_trials=tic;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Legacy pag and polymer distribution generation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% 4..2.2 Pag and polymer positions generation
        
        posPAG      = randPosGen(SYSTEM_LIMITS,RHO_PAG*prod(SYSTEM_SIZE));
        posPolymer  = randPosGen(SYSTEM_LIMITS,RHO_POLYMER*prod(SYSTEM_SIZE));
        
        % Initializing the pag info structure (pagdata)
        
        pagdata.posPAG              =   posPAG;
        pagdata.posPAG_removed      =   [];
        pagdata.acid_act_xyz_idx    =   []; % The index of the voxel
        pagdata.acid_act_e_xyz      =   []; % The position of event
        
        % Initializing the polymer info structure (polymdata)
        
        polymdata.posPolymer        =   posPolymer;
        polymdata.SE_act_xyz        =   [];
        polymdata.moleculeDensity   =   MOLECULAR_NUMBER_DENSITY;
        
        xyz_electron=[];

        totalacids_pertrial=0;        
        nacids_unsat_total=0;        
        
        %% 4..3 Loop through all incidence of electrons
        
        %%% The sequence of positions to be dosed
        dosingSequence     =   dosingSequenceHandle(nElectrons);
        energyScanArchive{trial_count}.dosingSequence = dosingSequence;
                
        %%% Initializing the archive element (delayed till this point to
        %%% accomodate stochastic dosing (different number of electron per
        %%% each trial
        energyScanArchive{trial_count}.incidences    =   cell([1 size(dosingSequence,2)]);
        
        for incidence   =   1:size(dosingSequence,2)
            %% 4..3.1 Iteration for dosed positions
            tst_incidence=tic;
            %fprintf(logfile_fid,'\n...Co-ordinate %d of %d\n',...
            %    incidence,size(dosingSequence,2));

            xval = dosingSequence(1,incidence);
            yval = dosingSequence(2,incidence);
            zval = dosingSequence(3,incidence);

            %fprintf('......Electron %d of %d\n',incidence,nElectrons);
            %fprintf(logfile_fid,'......Electron %d of %d\n',incidence,nElectrons);
            eventP   =   eventPrototype;
            eventP.xyz  =   [xval yval zval];

            eventP.theta_in = thetaGenHandle(1);
            theta_init  =   [theta_init eventP.theta_in];
            rng('shuffle');
            eventP.phi_in=2*pi*rand; 

            %%% Initialize the array that stores all events from this
            %%% electron
            eventdata={};

            % Again, registering the pre scattering pag 
            % and polymer loading
            posPAG_init     =   pagdata.posPAG;
            posPAG_rmv_init =   pagdata.posPAG_removed;
            posPolymer_init =   polymdata.posPolymer;

            %% 4..3.2The scattering engine itself
            [eventdata,pagdata,polymdata]=TrajectoryWrapper({eventP},scattdata,eventdata,xyzglobal,pagdata,polymdata,logfile_fid);

            % Display the location of acids created by this electron.
            %{
            if echoConfig.acid.perElectron
                fprintf('Positions of acids activted by this electron\n');
                disp(pagdata.posPAG_removed(:,...
                    size(posPAG_rmv_init,2)+1:...
                    size(pagdata.posPAG_removed,2))');
            end
            %}
            
            %% 4..3.3 Documenting the results
            
            energyScanArchive{trial_count}.incidences{incidence} = eventdata;

            % Record the time spent on each incidence
            tend_coordinates=toc(tst_incidence);
            %fprintf(logfile_fid,'...Took %.2f s\n',tend_coordinates);
            
        end
        
        %% 4..2.3 Post Trial Logging and counting
        
        %%% Do acid counting etc. right before next trial and energy values
        %%% are used:
        %%% The indices of the posisions where acid is generated
        acid_act_xyz_idx    =   pagdata.acid_act_xyz_idx;
        %Tranpose for compatibility with comparison code
        acid_xyz            =   pagdata.posPAG_removed';
        acid_act_e_xyz      =   pagdata.acid_act_e_xyz;
        %%% For consistency with legacy naming
        acid_fine_xyz       =   acid_act_e_xyz;
        
        %%% Acid counting for the trial
        nAcids_thruTrial(trial_count)       =   size(acid_xyz,1);
        meanAcids_thruTrial(trial_count)    =   size(acid_xyz,1)/size(dosingSequence,2);
        
        
        %%% Acid echo
        if echoConfig.acid.perTrial
            fprintf('Positions of acids activted in this trial (acid_xyz)\n');
            disp(acid_xyz)
            fprintf('Positions of acids activted in this trial (NaN counting)\n');
            disp(posPAG(:,isnan(pagdata.posPAG(1,:)))')
            fprintf('Positions of acids activted in this trial (array_indexing)\n');
            disp(posPAG(:,acid_act_xyz_idx)')
        end
        
 
        %%% Same thing for ions
        SE_act_xyz  =   polymdata.SE_act_xyz;
        %%% The acuulmulative arrays
        %ionXyz_TrAccu     =   [ionXyz_TrAccu; SE_act_xyz'];
        nIons_thruTrial(trial_count)        =   size(SE_act_xyz,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Putting things into the database
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% The positions of the acid and activation
        energyScanArchive{trial_count}.acid_xyz       = acid_xyz;
        energyScanArchive{trial_count}.activation_xyz = acid_act_e_xyz;
        % The array indices are not included as they provide no information
        % unless I includ the acid array as well
        energyScanArchive{trial_count}.ion_xyz        = SE_act_xyz;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    %% 4..1.3 Per Energy House Keeping
    meanAcids_thruE(E_count)    =   mean(nAcids_thruTrial);
    stdAcids_thruE(E_count)     =   std(nAcids_thruTrial);
    
    scanArchive{E_count}        =   energyScanArchive;
    
    %%% 
    %%% Legacy variable length array
    for trialIter = 1:nTrials
        posAcid_TrAccu      =   [posAcid_TrAccu;...
            energyScanArchive{trialIter}.acid_xyz];
        posAcidAct_TrAccu   =   [posAcidAct_TrAccu;...
            energyScanArchive{trialIter}.activation_xyz];
    end
    
    tend    = toc(tStart);    
    tEnd_E  = toc(tStart_E);
    fprintf(logfile_fid,'Total time elapsed = %.4f s = %.4f min = %.4f hours\n',tend,tend/60,tend/3600);
    fprintf(logfile_fid,'Of which %.4f s = %.4f min = %.4f hours was spent ib the last energy cycle\n',...
        tEnd_E,tEnd_E/60,tEnd_E/3600);
    
    fprintf('Total time elapsed = %.4f s = %.4f min = %.4f hours\n',tend,tend/60,tend/3600);
    fprintf('Of which %.4f s = %.4f min = %.4f hours was spent on the last energy cycle\n',...
        tEnd_E,tEnd_E/60,tEnd_E/3600);
    fprintf(logfile_fid,'%d; <Acids> = %.4f; sig_Acids = %.4f\n',nTrials,mean(nAcids_thruTrial),std(nAcids_thruTrial));
    fprintf(1,'%d trials; <Acids> = %.4f; sig_Acids = %.4f\n',nTrials,mean(nAcids_thruTrial),std(nAcids_thruTrial));    
    
    %%% Archive saved in each energy step in case of unexpected termination
    %{
    save(sprintf(strcat(outputBasePath,...
            strcat('ScanArchive_',num2str(eSweep(E_count))))),'energyScanArchive');
    %}
    % Conversion to a newer data base format
    parfor ii = 1:nTrials
        currDB(ii)    =   reorderArchive(energyScanArchive(ii))
    end
    databaseName = strcat('DB_',num2str(eSweep(E_count)),'eV');
    saveDatabase(currDB,outputBasePath,databaseName);
    % Exporting the archive
    archiveName = strcat('ScanArchive_',num2str(eSweep(E_count)),'eV');
    saveArchive(energyScanArchive,outputBasePath,archiveName);
    
    %% 4..1.4 Per Energy Trajectory Echo
    if echoConfig.traj3.perEnergy
        edHandle = figure(round(FIGURE_TRAJ_PER_ENERGY + eSweep(E_count)));
        edHandle.Name = strcat('Trajectory-',num2str(eSweep(E_count)),'eV');
        %%% Count the number of events
        totalEventCountAtEnergy = 0;
        grid on
        hold off
        qX = [];
        qY = [];
        qZ = [];
        qU = [];
        qV = [];
        qW = [];
        for trialIter   = 1:nTrials
            for inIterator = 1:size(scanArchive{E_count}{trialIter}.incidences,2)
                nEvents     =   length(...
                    scanArchive{E_count}{trialIter}.incidences{inIterator});
                for eventIterator = 1:nEvents
                    totalEventCountAtEnergy = totalEventCountAtEnergy +1;
                    this = scanArchive{E_count}{trialIter}.incidences{inIterator}{eventIterator};
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
        title({strcat('Acid and Trajectory visuallization at KE = ',num2str(eSweep(E_count)),' eV');...
            strcat(num2str(nTrials),' trials with',' ',...
            num2str(nElectrons),' electrons per trial on average')});
        legend({'Trajectories','Activations','Acids'});
        daspect([1 1 1]);
        drawnow;
    else
        %%% Still the number of events needs to be counted for the legacy
        %%% variables. Should have moved them to outside of the loop
        totalEventCountAtEnergy = 0;
        for trialIter   = 1:nTrials
            for inIterator = 1:size(scanArchive{E_count}{trialIter}.incidences,2)
                nEvents     =   length(...
                    scanArchive{E_count}{trialIter}.incidences{inIterator});
                for eventIterator = 1:nEvents
                    totalEventCountAtEnergy = totalEventCountAtEnergy +1;
                end
            end
        end
    end    
    
    %% 4..1.5 Per Energy Acid Distribution Echo
    if echoConfig.acidDist.active
        dispR = echoConfig.acidDist.range;
        dispS = echoConfig.acidDist.res;
        
        binEdges = { -dispR-dispS/2:dispS:dispR+dispS/2, -dispR-dispS/2:dispS:dispR+dispS/2 };
        axislabels = {'x (nm)','y(nm)','z(nm)'};
        
        esHandle    =   figure(round(FIGURE_ACID_DIST_PERENERGY + eSweep(E_count)));
        esHandle.Name = strcat('Summary-',num2str(eSweep(E_count)),'eV');
        clf; 
        
        annotation('textbox',[0.15 0.0025 .23 .05],'String','Linear Acid Density Maps')
        annotation('textbox',[0.7 0.9475 .21 .05],'String','Log Acid Density Maps')
        for ii = 1:3
            for jj = ii+1:3
                [freq,posC] = hist3(posAcid_TrAccu(:,[jj ii]),binEdges);
                sp = subplot(3,3,sub2ind([3 3],ii,jj));
                hold off
                im = imagesc([min(posC{2}) max(posC{2})],...
                    [min(posC{1}) max(posC{1})],...
                    freq);
                daspect([1 1 1])
                colormap(sp,'winter');
                hold on
                contour(posC{2},posC{1},...
                    freq,[0 0.05*max(max(freq))],'w')
                ylabel(axislabels{jj})
                xlabel(axislabels{ii})

                sp = subplot(3,3,sub2ind([3 3],jj,ii));
                hold off;
                imagesc([min(posC{2}) max(posC{2})],...
                    [min(posC{1}) max(posC{1})],...
                    log(freq));
                daspect([1 1 1])
                hold on
                colormap(sp,'default');
                contour(posC{2},posC{1},...
                    log(freq),3,'w')
                ylabel(axislabels{jj})
                xlabel(axislabels{ii})
            end
            
            subplot(3,3,sub2ind([3 3],ii,ii))
            hold off;
            histogram(posAcid_TrAccu(:,ii),binEdges{1},...
                'Normalization','pdf');
            hold on;
            histogram(posAcidAct_TrAccu(:,ii),binEdges{1},...
                'Normalization','pdf');
            xlabel(axislabels{ii})
            ylabel('Probability Density');
            if ii ==1
                legend({'Acid','Activation'},...
                    'Position',[0.2 0.9 0.1 0.025],...
                    'FontSize',6)
            end
        end 
    end        
    
end

%% 5 Post processing for backward compatibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The diagnostic script relies on certain arrays that I have abadoned. I
%%% will recreate them to avoid problems further down the road.
%%% === WARNING: THESE ARRAYS ONLY STORE DATA FROM THE LAST ENERGY ===
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_global    = zeros([1 totalEventCountAtEnergy]);
phi_global      = zeros([1 totalEventCountAtEnergy]);
stepLen_global  = zeros([1 totalEventCountAtEnergy]);
xyz_electron_global     = zeros([totalEventCountAtEnergy 3]);

jj = 0;
for trialIter   = 1:nTrials
    for inIterator = 1:size(scanArchive{E_count}{trialIter}.incidences,2)
        nEvents     =   length(...
            scanArchive{E_count}{trialIter}.incidences{inIterator});
        for eventIterator = 1:nEvents
            jj = jj + 1;
            this = scanArchive{E_count}{trialIter}.incidences{inIterator}{eventIterator};
            theta_global(jj)    =   this.theta;
            phi_global(jj)      =   this.phi;
            stepLen_global(jj)  =   this.pathlen;
            xyz_electron_global(jj,:) =   this.xyz;
        end
    end
end

%% 6.1 Additional Graphs
aveHandle = figure(FIGURE_ACID_VS_ENERGY);
aveHandle.Name = 'Acid_IncidentEnergy';
plot(eSweep,meanAcids_thruE,'-o','linewidth',3.0,'markersize',12);
errorbar(eSweep,meanAcids_thruE,stdAcids_thruE/sqrt(nTrials),...
    '--o','linewidth',2.0,'markersize',5);
xlabel('E (eV)');
ylabel('Mean number of acids');
set(gca,'fontsize',20,'linewidth',2.0);

aveHandle = figure(FIGURE_ACID_VS_ENERGY+1);
aveHandle.Name = 'Acid_IncidentEnergy_perEnergy';
plot(eSweep,meanAcids_thruE./eSweep,...
    '-o','linewidth',3.0,'markersize',12);
errorbar(eSweep,meanAcids_thruE./eSweep,...
    stdAcids_thruE./eSweep/sqrt(nTrials),...
    '--o','linewidth',2.0,'markersize',5);
xlabel('E (eV)');
ylabel('Mean number of acids per incident energy');
set(gca,'fontsize',20,'linewidth',2.0);


%% 6.2 Saving Graphs
if ~exist(strcat(outputBasePath, '\graphics'),'file')
       mkdir(strcat(outputBasePath, '\graphics'));
end
graphicFolderName = strcat(outputBasePath, '\graphics');   
figList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(figList)
  figHandle = figList(iFig);
  figName   = get(figHandle, 'Name');
  %savefig(figHandle, fullfile(graphicFolderName, figName, '.fig'));
  if ~isempty(num2str(figHandle.Number))>0
    saveas(figHandle,fullfile(graphicFolderName,...
        strcat('figure',num2str(figHandle.Number),...
        figName,'.pdf')));
  end
end

%% 6.2 Final workspace export
save(sprintf(strcat(outputBasePath,...
            'WorkSpace')));
diary off