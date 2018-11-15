%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What it does: start a trajectory with the secondary from the input event
%%%     scatter -> propagate -> scatter -> ..... and so on
%%% This is a toy model for calibration. The electron is assumed to have a
%%% fixed IMFP (the low energy IMFP) and each collision will cost a fixed
%%% amount of energy loss of 1 eV . By tunning the initial energy one
%%% tune the spread. The process is terminated by the stone wall kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [events,pagdata,polymdata]=TrajectoryFollowerRandomWalk(event,scattdata,scatt_Elim,xyzglobal,pagdata,polymdata,~)
%% Initializations
global illustration scattVector thetaLog;

    % Energy thresholds for various events
    pag_Eamin       =   5;
    E_ionize_min    =   12;
    E_inel_thr      =   scattdata.E_inel_thr;


    Eold            =   event.Ese;
    xold            =   event.xyz(1);
    yold            =   event.xyz(2);
    zold            =   event.xyz(3);
    theta_old       =   event.theta_in; % initial condition for theta in z
    phi_old         =   event.phi_in; % initial condition for the polar angle

    Enew=Eold;
    count=1;
    pagact_prob=0.1;

    %%%% get the PAG data
    posPAG              =   pagdata.posPAG;
    %position of removed pags
    posPAG_removed      =   pagdata.posPAG_removed;
    %%% The indices for querying the acid position array
    acid_act_xyz_idx    =   pagdata.acid_act_xyz_idx;
    %%% The location of the activation event
    acid_act_e_xyz      =   pagdata.acid_act_e_xyz;

    %%% polymer distribution matrix
    posPolymer      =   polymdata.posPolymer;
    SE_act_xyz      =   polymdata.SE_act_xyz;

    %%% Reaction Radius
    pag_rcnrad      =   event.pag.rcnrad;
    %%% Molecule density per cm^3
    moleculeDensity      =   polymdata.moleculeDensity;
    
    %% The steps loop
    while Eold>scatt_Elim
        %% 2. Calculating the IMFP
        %%% It's a bit tricky here. So the engine has the option to give
        %%% partial only IMFP, only IMFP is calculated. If not the engine
        %%% will be run from top to bottom
        
        %% 2.1.1 Random walk
        imfp_randomWalk     =   event.lowEimfp;
        
        %% 2.1.3 Stone Wall type low energy cutoff
        stoneWallResults = scattEngineStoneWall(Eold,scattdata);
        imfp_stoneWall      =   stoneWallResults.imfp;
        theta_stoneWall     =   stoneWallResults.theta;
        phi_stoneWall       =   stoneWallResults.phi;
        eLoss_stoneWall     =   stoneWallResults.eLoss;
        rxnRadius_stoneWall =   stoneWallResults.rxnR;
        
        %% 2.2 Caculate the total IMFP
        imfp    =   1/(1/imfp_randomWalk+1/imfp_stoneWall);            
        
        %% 3. Electron propagation, chronologically before the scattering
        rnew=exprnd(imfp); % exponential distribution        

        znew=zold+rnew*cos(theta_old);
        xnew=xold+rnew*sin(theta_old)*cos(phi_old);
        ynew=yold+rnew*sin(theta_old)*sin(phi_old);

        %%%% calculate the path length
        pathlen=sqrt((xnew-xold)^2+(ynew-yold)^2+(znew-zold)^2);
        
        %% 4. Scattering Logic
        
        %% 4.1  Determin which collision mechanism is at play
        
        %%% The probability (or should I say odds)of each event is
        %%% proportional to their inverse IMFP. The original code was right
        %%% but I'm adding a more genral framework in case additional
        %%% scattering mechanisms come along
        
        interactionCandidate  =   {'RandomWalk','StoneWall'};
        scattInvIMFP    =   [(1/imfp_randomWalk) (1/imfp_stoneWall)];
        
        scattType = weightedCategoricalRandGen(interactionCandidate,scattInvIMFP);
        
        switch(scattType)
            case 'RandomWalk'
                act         =   'RandomWalk';
                Eloss_val   =   1;
                theta       =   acos(2*rand-1);
                phi         =   2*rand*pi;
                rxnRadius   =   0;
            case 'StoneWall'
                act         =   'StoneWall';
                Eloss_val   =   eLoss_stoneWall;
                theta       =   theta_stoneWall;
                phi         =   phi_stoneWall;
                rxnRadius   =   rxnRadius_stoneWall;
        end        
        
        %% 4.2  Transfomration back into the resist frame
        % The angles in the scattering results are relative to the
        % pre-event travelling direction. Thus a coordinate transformation
        % is needed
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Debug: To select a known form of scattering
        %%%================================================================
            scatteringDebug = 'none';
        %%% (options are 'forward', 'orthogonal', 'random' and 
        %%% 'none'(which is equal to not intervening)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch (scatteringDebug)
            case 'forward'
                theta   =   0;
                phi     =   rand*2*pi;
            case 'orthogonal'
                theta   =   pi/2;
                phi     =   rand*2*pi;
            case 'random'
                theta   =   acos(2*rand-1);
                phi     =   rand*2*pi;
            otherwise
        end                
        
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Tracking changes in direction
        scattVector =   [scattVector newUnitVec];
        thetaLog    =   [thetaLog   theta_old];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% 4.3  Change in direction and updating the angle iterator
        vec2D_mag   =   sqrt(newUnitVec(1).^2+newUnitVec(2).^2);
        theta_new   =   atan2(vec2D_mag,newUnitVec(3));
        phi_new     =   atan2(newUnitVec(2),newUnitVec(1));
                
        %% 5. Energy Deposition [PAG activations? SE-gen?]
        %%% In react-propagate scenario, scattering takes place at the "old
        %%% coordinates"
        %%% In propagate-react scenario, scattering takes place at the "new
        %%% coordinates"
        
        xEvent     =   xnew;
        yEvent     =   ynew;
        zEvent     =   znew;
        
        [pagidx,npags,polymidx,npolyms]=...
            pag_polym_query([xEvent yEvent zEvent],posPAG,posPolymer,rxnRadius);
        
        pag_ratio=npags/(npags+npolyms);
        npags_removed=0;
        
        %% 5.2 The possibilities
        %{
        if (rand<pag_ratio...%
                || strcmp(act,'6eVRes')... %
                || strcmp(act,'LowEnergy-Acid'))...% if activating a PAG
                && ~(strcmp(act,'vibr')&& Eloss_val<pag_Eamin)  % Vibrational excitations treated differently 
         %% 5.1.1 Acid generation (by volume ratio and 6eV resonance)
            Ese=0;
            nSE=0;
            nion=0;
            nacid=0;
            nacid_unsat=1;            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%  Satuation test 
            satuationTest = 0;
            if satuationTest
                nacid=1; 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            if npags>0
                num2Remove=1; % how many PAGS to remove [could be Eloss/Eact in the future!]
                npags_removed=npags_removed+1;
                nacid=1;
                %%% Determine which pag to remove 
                [posPAG,posPAG_removed,acid_act_xyz_idx,acid_act_e_xyz]...
                    =acidActivation(num2Remove,posPAG,posPAG_removed,...
                    pagidx,[xEvent, yEvent, zEvent],...
                    acid_act_xyz_idx,acid_act_e_xyz);                
                
                %%% Determine if this event spawns a secondary electron
                Ese=Eloss_val-pag_Eamin;    
                Ese(Ese<0)=0;
                if Ese>0
                    nSE=1;
                    nion=1;
                    SE_act_xyz=[SE_act_xyz posPAG_removed];
                end

                act='acid';                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%     Debug code to figure out why the fine position arrays is
                %%%     not as long as the pixelated one
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %{
                if length(acid_act_xyz_idx)~=size(acid_act_xyz,1)
                    debugObject = {};
                    debugObject.condition   =   'The two acid vectors have different lengths';
                    debugObject.level       =   'trajcalc';
                    debugObject.coarseArrayLength   = length(acid_act_xyz_idx);
                    debugObject.fineArrayLength     = size(acid_act_xyz,1);
                    debugOutput{size(debugOutput,2)+1} = debugObject;
                end                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %}
            else
    %             Eloss_val=0; % if no pag available, no energy loss.
    %             nacid=0;
                E_benzene=6.5;
                Ese=Eloss_val-E_benzene;    Ese(Ese<0)=0;
                act='acid-none-Polym-noSE';
                if Ese>0 && ~isempty(polymidx)
                    %%%Signaling the engine to initiate a secondary
                    nSE=1;
                    nion=1;
                    num2add=nion;
                    add_idx=randi([1 length(polymidx)],num2add);
                    SE_act_xyz=[SE_act_xyz posPolymer(:,polymidx(add_idx))];
                    act='acid-none-Polym-SE';
                end
            end
        else
            nacid=0;
            nacid_unsat=0;
            nion=0;
            nSE=0;
            Ese=0;
            switch act
                case 'SE'
                %% 5.1.2 Secondary electron generation
                    gen_SE=1; % if 1, it will generate secondaries, if 0 then no SE is generated
                    if npolyms>0 % if there was no polymer molecule
                        Ese=(Eloss_val>E_ionize_min).*(Eloss_val-E_ionize_min);
                        nSE=double(Eloss_val>E_ionize_min); % set to 1 if above is set up for SE-gen instead of creating excited states
                        nion=nSE;
                        num2add=nion;
                        if num2add>0 && ~isempty(polymidx)
                            add_idx=randi([1 length(polymidx)],num2add);
%                             polym_img(polymidx(add_idx))=polym_img(polymidx(add_idx))+1;
                            SE_act_xyz=[SE_act_xyz posPolymer(:,polymidx(add_idx))];
                        end                        
                    else
                        Eloss_val=0; % set to 0 if there was no event
                    end
                case '6eVRes'   %%% !!!! This branch is not reachable
                    act='6eVRes-Polym';
                case 'vibr'     %%% !!! This 
                    %nacid=1;% temporary, added on 5/27/2017 to test saturation 
                    nacid_unsat=1;
                    
                    %%% Acid generation is mendatory. That looks dubious.
                    %%% -JHM
                    if npags>999
                        
          %             Eloss_val=min([6.8 max([Eloss_val pag_Eamin])]); % modify this value, as its a 6.8 eV event instead.
                        num2Remove=1; % how many PAGS to remove [could be Eloss/Eact in the future!]
                        npags_removed=npags_removed+1;
                        nacid=1;
                        
                        [posPAG,posPAG_removed,acid_act_xyz_idx,acid_act_e_xyz]...
                            =acidActivation(num2Remove,posPAG,posPAG_removed,...
                            pagidx,[xEvent, yEvent, zEvent],...
                            acid_act_xyz_idx,acid_act_e_xyz);                        
                        
                        Eloss_val=pag_Eamin;
                        Ese=Eloss_val-pag_Eamin;    Ese(Ese<0)=0;
                        if Ese>0
                            nSE=1;
                            nion=1;
                            SE_act_xyz=[SE_act_xyz posPolymer(:,pagidx(remove_idx))];
                        end
                        acid_act_xyz_idx=[acid_act_xyz_idx pagidx(remove_idx)];
                        act='vibr-acid';
                    else
                        act='vibr-polym';
                    end
                    
                    %%% Conventional SE logic. Tempory -JHM
                    Ese     =   (Eloss_val>E_ionize_min).*(Eloss_val-E_ionize_min);
                    nSE     =   double(Eloss_val>E_ionize_min);
                    
                    Ese(Ese<0)=0;
                    if Ese>0
                        nSE=1;
                        nion=1;
                        SE_act_xyz=[SE_act_xyz posPolymer(:,pagidx(remove_idx))];
                    end
                case 'StoneWall'
                    num2Remove  =   1;
                    [posPAG,posPAG_removed,acid_act_xyz_idx,acid_act_e_xyz]...
                        =acidActivation(num2Remove,posPAG,posPAG_removed,...
                        pagidx,[xEvent, yEvent, zEvent],...
                        acid_act_xyz_idx,acid_act_e_xyz);
            end
        end
        %}
        
        if strcmp(act, 'StoneWall')
            num2Remove  =   1;
                    [posPAG,posPAG_removed,acid_act_xyz_idx,acid_act_e_xyz]...
                        =acidActivation(num2Remove,posPAG,posPAG_removed,...
                        pagidx,[xEvent, yEvent, zEvent],...
                        acid_act_xyz_idx,acid_act_e_xyz);
            nacid   =   1;            
        else
            nacid   =   0;
        end
        
        nacid_unsat = -1;
        
        Enew    =   Eold-Eloss_val;
        nSE     =   0;
        Ese     =   0;
        
        %% 6. Documenting the event
        xyzglobal.x=[xyzglobal.x xnew];
        xyzglobal.y=[xyzglobal.y ynew];
        xyzglobal.z=[xyzglobal.z znew];
        
        %%%% post energy-loss variables:
        events{count}.xyz_init  =   [xold yold zold];
        events{count}.xyz       =   [xnew ynew znew];
        events{count}.xyzglobal =   xyzglobal;
        events{count}.pathlen   =   pathlen;
        events{count}.Ein       =   Eold;
        events{count}.Eout      =   Enew;
        events{count}.Eloss     =   Eloss_val;
        events{count}.act       =   act;
        events{count}.imfp      =   imfp;
        events{count}.rnew      =   rnew;   
        %%% parameters in scattering frame and the scattering
        events{count}.theta     =   theta;
        events{count}.phi       =   phi;
        events{count}.scattType =   scattType;
        %%% parameters in resist frame
        events{count}.theta_in  =   theta_old;
        events{count}.theta_out =   theta_new;
        events{count}.phi_in    =   phi_old;
        events{count}.phi_out   =   phi_new;

        events{count}.Ese       =   Ese;
        events{count}.nacid     =   nacid;
        events{count}.nacid_unsat   =   nacid_unsat;
        events{count}.nSE       =   nSE;

        events{count}.pag.rcnrad=pag_rcnrad;
        
        %% 7. Updating the "old" (initial) variable values
        Eold    =   Enew;

        xold        =   xnew;
        yold        =   ynew;
        zold        =   znew;
        theta_old   =   theta_new;
        phi_old     =   phi_new;

        count=count+1;
        
    end
    %% 8. Update the resist composition
    pagdata.posPAG              =   posPAG;
    pagdata.posPAG_removed      =   posPAG_removed;
    pagdata.acid_act_xyz_idx    =   acid_act_xyz_idx;
    pagdata.acid_act_e_xyz      =   acid_act_e_xyz;

    pagdata.posPolymer          =   posPolymer;
    polymdata.SE_act_xyz        =   SE_act_xyz;
end
