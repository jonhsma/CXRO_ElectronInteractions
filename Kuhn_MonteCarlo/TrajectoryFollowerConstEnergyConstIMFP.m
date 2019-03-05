% This function followes a trajectory in the Monte Carlo electron trajectory 
% model with no energy loss and energy deposition. 
% A few key differeces from the one used in the discrete energy loss sim:
% 1. The length of the trajectory specified
% 2. Everything is stored in arrays to speed things up


function segments = TrajectoryFollowerConstEnergyConstIMFP(initConfig,...
    scattdata,... The scattering data
    THFmoleculeDensity,... Tempory
    structureArrayAllocator... To make sure that I can concatenate the arrays
    )
%% 1 Initializations

    enenergy        =   initConfig.energy;
    x_old           =   initConfig.xyz(1);
    y_old           =   initConfig.xyz(2);
    z_old           =   initConfig.xyz(3);
    theta_old       =   initConfig.theta_in; % initial condition for theta in z
    phi_old         =   initConfig.phi_in; % initial condition for the polar angle
    
    nSteps          =   initConfig.nSteps;
    
    %The fixed imfp
    FIXED_IMFP = 0.5;


    %%% Molecule density per cm^3
    %%% This is a tempory solution
    moleculeDensity      =   THFmoleculeDensity;
    
    % Allocate memory for speed
    segments = structureArrayAllocator(nSteps);
    
    % Some adaptations
    event.lowEimfp=3.67;
    
    %% The steps loop
    for ii = 1:1:nSteps
        %% 2. Calculating the IMFP
        %%% It's a bit tricky here. So the engine has the option to give
        %%% partial only IMFP, only IMFP is calculated. If not the engine
        %%% will be run from top to bottom        
        % 2.1.1 The Optical component that always gets calculated
        imfp_opt = genMFP_OptData(event,scattdata.optical,enenergy);
        % 2.1.2 Vibrational calculations
        imfp_vibr = genMFP_Vibr(scattdata.vibr,enenergy,moleculeDensity);     
        % 2.1.3 Stone Wall type low energy cutoff
        imfp_stoneWall      =   genMFP_StoneWall(enenergy,scattdata.stoneWall);        
        % 2.2 Caculate the total IMFP
        imfp    =   FIXED_IMFP;                    
        %% 3. Electron propagation, chronologically before the scattering
        deltaR=exprnd(imfp); % exponential distribution     

        z_new=z_old+deltaR*cos(theta_old);
        x_new=x_old+deltaR*sin(theta_old)*cos(phi_old);
        y_new=y_old+deltaR*sin(theta_old)*sin(phi_old);

        %%%% calculate the path length, this is an redundancy
        pathlen=sqrt((x_new-x_old)^2+(y_new-y_old)^2+(z_new-z_old)^2);
          
        %% 4. Scattering Logic
        
        %%% The probability (or should I say odds)of each event is
        %%% proportional to their inverse IMFP. The original code was right
        %%% but I'm adding a more genral framework in case additional
        %%% scattering mechanisms come along
        
        interactionCandidate  =   {'Optical','Vibrational','StoneWall'};
        scattInvIMFP    =   [(1/imfp_opt) (1/imfp_vibr) (1/imfp_stoneWall)];
        
        scattType = weightedCategoricalRandGen(interactionCandidate,scattInvIMFP);
        
        switch(scattType)
            case 'Vibrational'
                vibrScattCmplx     =   genrandEloss_Vibr(scattdata.vibr,enenergy);
                Eloss_val   =   vibrScattCmplx.Eloss;
                theta       =   vibrScattCmplx.theta;
                phi         =   vibrScattCmplx.phi;
                act         =   'vibr';
            case 'Optical'
                controlparm.onlyimfp=0;
                if isfield(scattdata.optical,'inel_dcsdata')
                    Elossrand_opt=genrandEloss_OptData_JHM(scattdata.optical,enenergy,scattdata.optical.inel_dcsdata,controlparm);
                else
                    Elossrand_opt=genrandEloss_OptData_JHM(scattdata.optical,enenergy,controlparm);
                end
                Eloss_opt   =   Elossrand_opt.Eloss; % the optical energy loss
                theta_opt   =   Elossrand_opt.theta; % the optical theta
                phi_opt     =   Elossrand_opt.phi; % the optical phi

                Eloss_val=Eloss_opt;
                theta = theta_opt;
                phi=phi_opt;
                act = 'Optical';
                
            case 'LowEnergy' % Low energy random walk, Currently unaccessible 
                act='LowEnergy';
                theta = acos(2*rand-1);
                phi=2*pi*rand;  
            case 'StoneWall'
                act         =   'StoneWall';
                stoneWallResults = scattEngineStoneWall(enenergy,scattdata.stoneWall);
                Eloss_val   =   stoneWallResults.eLoss;
                theta       =   stoneWallResults.theta;
                phi         =   stoneWallResults.phi;
                rxnRadius   =   stoneWallResults.rxnR;
        end        
        
        %% 4.2  Transfomration back into the resist frame
        % The angles in the scattering results are relative to the
        % pre-event travelling direction. Thus a coordinate transformation
        % is needed
        
        %%% Transforming the results of the scattering simulation back into
        %%% the lab coordinates
        
        cT_ini = cos(theta_old);
        sT_ini = sin(theta_old);
        
        cP_ini = cos(phi_old);
        sP_ini = sin(phi_old);
        
        rToResist_theta =...
            [   cT_ini,  0,     sT_ini;...
                0,       1,     0;...
                -sT_ini,       0,    cT_ini];
            
        rToResist_phi =...
            [   cP_ini, -sP_ini     0;...
                sP_ini, cP_ini,     0;...
                0,          0,      1];
            
        collisionFrameVector =...
            [sin(theta)*cos(phi);...
            sin(theta)*sin(phi);...
            cos(theta)];
        
        newUnitVec = rToResist_phi*rToResist_theta*collisionFrameVector; 
        
        %% 4.3  Change in direction and updating the angle iterator
        vec2D_mag   =   sqrt(newUnitVec(1).^2+newUnitVec(2).^2);
        theta_new   =   atan2(vec2D_mag,newUnitVec(3));
        phi_new     =   atan2(newUnitVec(2),newUnitVec(1));        
        %% 6. Documenting the event
        
        % Post propagation values
        segments(ii).xyz_init  =   [x_old y_old z_old]';
        segments(ii).xyz_final =   [x_new y_new z_new]';
        segments(ii).xyz_delta =   [x_new y_new z_new]'-[x_old y_old z_old]';
        % Position vectors are transposed for easy data processing
        segments(ii).pathlen   =   pathlen;
        segments(ii).Eloss     =   Eloss_val;
        segments(ii).imfp      =   imfp;
        segments(ii).deltaR    =   deltaR;   
        %%% parameters in scattering frame or San Francisco
        segments(ii).theta_SF       =   theta;
        segments(ii).phi_SF         =   phi;
        segments(ii).scattType      =   scattType;
        segments(ii).act            =   act;
        %%% parameters in lab frame
        segments(ii).theta_init  =   theta_old;
        segments(ii).theta_final =   theta_new;
        segments(ii).phi_init    =   phi_old;
        segments(ii).phi_final   =   phi_new;
        
        %% 7. Updating the "old" (initial) variable values
        x_old        =   x_new;
        y_old        =   y_new;
        z_old        =   z_new;
        theta_old   =   theta_new;
        phi_old     =   phi_new;        
    end
end