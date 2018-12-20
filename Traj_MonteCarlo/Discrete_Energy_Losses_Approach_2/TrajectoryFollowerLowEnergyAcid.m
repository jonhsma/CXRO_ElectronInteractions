%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What it does: start a trajectory with the secondary from the input event
%%%     propagate -> scatter -> scatter ->..... and so on
%%% This version doesn't allow for acid generation from plasmon scattering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [events,pagdata,polymdata]=TrajectoryFollowerLowEnergyAcid(event,scattdata,scatt_Elim,xyzglobal,pagdata,polymdata,~)
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
        %% 1. Determine the active scattering mechanisms        
        if strcmp(scattdata.vibr.datasrc,'Khakoo')==1
            vibr_ics    =   scattdata.vibr.ics;
            vibr_ics(isnan(vibr_ics))...
                        = 0;
            vibrdata_E  =   vibr_ics(:,1);
            if Eold     >   max(vibrdata_E) % no vibrational data here, so go with inelasic
                scattType='HE'; % scattering is of "high-energy" type
                if illustration
                    fprintf('scattType = HE; count: %d \n',count);
                end
            else
                scattType='LE'; % consider "low-energy" as well
                if illustration
                    fprintf('scattType = LE; count: %d \n',count);
                end
            end
        end
        if strcmp(scattdata.vibr.datasrc,'Frohlich')==1
            scattType='HE'; % since its analytical, can evaluate it at an arbitrary energy
        end
        if ~exist('scattType')
            error('......TrajCalc2: The vibrational scattering data could not be read. Check the ''scattdata.vibr'' data structure');
        end        
        
        %% 2. Calculating the IMFP
        %%% It's a bit tricky here. So the engine has the option to give
        %%% partial only IMFP, only IMFP is calculated. If not the engine
        %%% will be run from top to bottom
        
        %% 2.1.1 The Optical component that always gets calculated
        controlparm.onlyimfp=1;
        if isfield(scattdata.optical,'inel_dcsdata')
            Elossrand_opt=genrandEloss_OptData_JHM(scattdata.optical,Eold,scattdata.optical.inel_dcsdata,controlparm);
        else
            Elossrand_opt=genrandEloss_OptData_JHM(scattdata.optical,Eold,controlparm);
        end
        if Eold>event.lowEthr
            controlparm.onlyimfp=0;
            imfp_opt=Elossrand_opt.imfp; % imfp in nm; the above lines only calculate imfp due to the controlparms.onlyimfp line abvove.
            if isnan(imfp_opt)
                imfp_opt=Inf;
            end
        else
            imfp_opt=event.lowEimfp;
        end
        
        %% 2.1.2 Vibrational calculations
        if strcmp(scattType,'LE')==1 % Khakoo data from above if-else
            if strcmp(scattdata.vibr.datasrc,'Khakoo')==1
                 vibrScattCmplx     =   genrandEloss_Vibr(scattdata.vibr,Eold);
                 eLoss_vibr         =   vibrScattCmplx.Eloss;
                 invImfp_invCM      =   vibrScattCmplx.ics*moleculeDensity;
                 imfp_vibr          =   1/invImfp_invCM*1e7; % imfp in nm
                 theta_vibr         =   vibrScattCmplx.theta;
                 phi_vibr           =   vibrScattCmplx.phi;
            else % Frohlich data from above if-else
                %{
                if strcmp(scattdata.vibr.datasrc,'Frohlich')==1
                    if Eold<=Inf % just temporary to speed things up [Inf if you want to execute whats in the if-block]
                        eps0=scattdata.vibr.eps0;
                        epsinf=scattdata.vibr.epsinf;
                        E_optphonon=scattdata.vibr.hbarw;
                        imfp_vibr1=[];
                        for Eph_count=1:length(E_optphonon)
                            [imfp_vibr,imfp_vibr_creation,imfp_vibr_annihilation,theta_vibr]=scattdata.vibr.imfp_func(Eold,eps0,epsinf,E_optphonon);
                            [imfp_vibr1(Eph_count),imfp_vibr_creation,imfp_vibr_annihilation,theta_vibr]=scattdata.vibr.imfp_func(Eold,eps0,epsinf,E_optphonon(Eph_count));
                            [imfp_vibr2,imfp_vibr_creation,imfp_vibr_annihilation,theta_vibr]=scattdata.vibr.imfp_func(Eold,eps0,epsinf,E_optphonon(Eph_count));
                        end
                        imfp_vibr=1/(sum(1./imfp_vibr1));
                        imfp_vibr=1/(1/imfp_vibr1 + 1/imfp_vibr2);
                        if rand<(1/imfp_vibr1)/(1/imfp_vibr1 + 1/imfp_vibr2)
                            Eloss_val=E_optphonon(1);
                        end
                        Eloss_vibr=E_optphonon(1);
                        Eloss_vec=scattdata.vibr.E;
                        imfp_vec=scattdata.vibr.imfp;
                        imfp_vibr=interp1(Eloss_vec,imfp_vec,Eold,'linear','extrap');
                        theta_vibr=rand*2*pi;
                        rng('shuffle');
                        phi_vibr=-pi+2*pi*rand;
                    else
                        imfp_vibr=Inf;
                        phi_vibr=0;
                        theta_vibr=0;
                    end
                end
                %}
                fprintf('Somehow it gets to the Frolich branch. Exiting\n');
            end
        else
            imfp_vibr=Inf;  % by default, no vibrational scattering
            phi_vibr=0;     % temporary, never gets used if the below is commented, so long as imfp_vibr=Inf as above
            theta_vibr=0;   % temporary, never gets used if the below is commented, so long as imfp_vibr=Inf as above
        end        
        if imfp_vibr==Inf && imfp_opt==Inf
%             fprintf(logfile_fid,'Warning in trajCalc3: IMFP (Vibr) = IMFP (Optical) = Inf\n');
        end         
        
        %% 2.1.3 Stone Wall type low energy cutoff
        stoneWallResults = scattEngineStoneWall(Eold,scattdata);
        imfp_stoneWall      =   stoneWallResults.imfp;
        theta_stoneWall     =   stoneWallResults.theta;
        phi_stoneWall       =   stoneWallResults.phi;
        eLoss_stoneWall     =   stoneWallResults.eLoss;
        rxnRadius_stoneWall =   stoneWallResults.rxnR;
        
        %% 2.2 Caculate the total IMFP
        imfp    =   1/(1/imfp_opt+1/imfp_vibr+1/imfp_stoneWall);            
        
        %% 3. Electron propagation, chronologically before the scattering
        rnew=exprnd(imfp); % exponential distribution        

        znew=zold+rnew*cos(theta_old);
        xnew=xold+rnew*sin(theta_old)*cos(phi_old);
        ynew=yold+rnew*sin(theta_old)*sin(phi_old);

        %%%% calculate the path length
        pathlen=sqrt((xnew-xold)^2+(ynew-yold)^2+(znew-zold)^2);
        
        %%% Check if the electron is still inside the resist
        %%% If not, count that as as escape
        limits  =   scattdata.SYSTEM_LIMITS;
        if ~(xnew <= limits(1,2)&& xnew >= limits(1,1) &&...
            ynew <= limits(2,2)&& ynew >= limits(2,1) &&...
            znew <= limits(3,2)&& znew >= limits(3,1))
        
            act = 'escape';
            
            %%%% post energy-loss variables:
            events{count}.xyz_init  =   [xold yold zold];
            events{count}.xyz       =   [xnew ynew znew];
            events{count}.xyzglobal =   xyzglobal;
            events{count}.pathlen   =   pathlen;
            events{count}.Ein       =   Eold;
            events{count}.Eout      =   Enew;
            events{count}.Eloss     =   0;
            events{count}.act       =   act;
            events{count}.imfp      =   imfp;
            events{count}.rnew      =   rnew;   
            %%% parameters in scattering frame and the scattering
            events{count}.theta     =   0;
            events{count}.phi       =   0;
            events{count}.scattType =   'escape';
            %%% parameters in resist frame
            events{count}.theta_in  =   theta_old;
            events{count}.theta_out =   theta_old;
            events{count}.phi_in    =   phi_old;
            events{count}.phi_out   =   phi_old;

            events{count}.Ese       =   0;
            events{count}.nacid     =   0;
            events{count}.nacid_unsat   =   0;
            events{count}.nSE       =   0;

            events{count}.pag.rcnrad=pag_rcnrad;
            
            break;
        end
        
        %% 4. Scattering Logic
        
        %% 4.1  Determin which collision mechanism is at play
        
        %%% The probability (or should I say odds)of each event is
        %%% proportional to their inverse IMFP. The original code was right
        %%% but I'm adding a more genral framework in case additional
        %%% scattering mechanisms come along
        
        interactionCandidate  =   {'Optical','Vibrational','StoneWall'};
        scattInvIMFP    =   [(1/imfp_opt) (1/imfp_vibr) (1/imfp_stoneWall)];
        
        scattType = weightedCategoricalRandGen(interactionCandidate,scattInvIMFP);
        
        switch(scattType)
            case 'Vibrational'
                %% 4.1.1 Vibrational Scattering
                %%% Engine has been ran and the results are in
                Eloss_val   =   eLoss_vibr;
                theta       =   theta_vibr;
                phi         =   phi_vibr;
                act         =   'vibr';
                rxnRadius = pag_rcnrad;
            case 'Optical'
                %% 4.1.1 Mermin based collision event simulation
                controlparm.onlyimfp=0;
                if isfield(scattdata.optical,'inel_dcsdata')
                    Elossrand_opt=genrandEloss_OptData_JHM(scattdata.optical,Eold,scattdata.optical.inel_dcsdata,controlparm);
                else
                    Elossrand_opt=genrandEloss_OptData_JHM(scattdata.optical,Eold,controlparm);
                end
                Eloss_opt   =   Elossrand_opt.Eloss; % the optical energy loss
                theta_opt   =   Elossrand_opt.theta; % the optical theta
                phi_opt     =   Elossrand_opt.phi; % the optical phi

                Eloss_val=Eloss_opt;
                theta = theta_opt;
                phi=phi_opt;
                
                %%% Set the reaction radius
                rxnRadius = pag_rcnrad;

                %%% Determine action by energy loss
                if Eloss_val==0
                    act='none'; % if Sp=0, no energy loss.
                elseif Eloss_val<pag_Eamin
                    act='Eloss<pagEaMin';
                elseif Eloss_val<E_ionize_min
                    act='6eVRes';
                else
                    act='SE';                                                    
                end
            case 'LowEnergy' % Currently unaccessible 
                %% 4.1.3 Low energy random walk (Currently inaccessible)
                act='LowEnergy';
                theta = acos(2*rand-1);
                phi=2*pi*rand;  
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

        %% 5.1 The possibilities
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
            
            if npags>0
                num2Remove=1; % how many PAGS to remove [could be Eloss/Eact in the future!]
                npags_removed=+1;
                nacid=1;
                %%% Determine which pag to remove 
                [posPAG,posPAG_removed,acid_act_xyz_idx,acid_act_e_xyz]...
                    =acidActivation(num2Remove,posPAG,posPAG_removed,...
                    pagidx,[xEvent, yEvent, zEvent],...
                    acid_act_xyz_idx,acid_act_e_xyz);                
                
                %%% Determine if this event spawns a secondary electron
                Ese         =   Eloss_val-pag_Eamin;    
                Ese(Ese<0)  =   0;
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
                E_benzene=6.5;
                Ese     =   Eloss_val-E_benzene;    Ese(Ese<0)=0;
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
        %}
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
                case 'vibr'   
                    act='vibr-polym';                    
                    
                    %%% Conventional SE logic. Tempory -JHM
                    %%% Should be removed as energy loss is always smaller
                    %%% than 1 eV
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
        %end

        Enew=Eold-Eloss_val;

        %fprintf(logfile_fid,'............trajFollower: [npags,len(pagidx),npags_removed,npolyms,pagratio,nacid,nion,nSE,Eloss,Ese,type] = [%d,%d,%d,%d,%.4f,%d,%d,%d,%.3f,%.3f,%s]\n',npags,length(pagidx),npags_removed,npolyms,pag_ratio,nacid,nion,nSE,Eloss_val,Ese,act);        
        
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
