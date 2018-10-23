%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What it does: start a trajectory with the secondary from the input event
%%%     scatter -> propagate -> scatter -> ..... and so on
%%% This is a toy model in an attempt to dix the anisotropy problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [events,pagdata,polymdata]=trajcalc3_JHSM_angleAdd(event,scattdata,scatt_Elim,xyzglobal,pagdata,polymdata,logfile_fid)
global illustration scattVector thetaLog;

    % Energy thresholds for various events
    pag_Eamin=5;
    E_ionize_min=12;
    E_inel_thr=scattdata.E_inel_thr;


    Eold=event.Ese;
    xold=event.xyz(1);
    yold=event.xyz(2);
    zold=event.xyz(3);
    theta_old=event.theta_in; % initial condition for theta in z
    phi_old=event.phi_in; % initial condition for the polar angle

    Enew=Eold;
    count=1;
    pagact_prob=0.1;

    %%%% get the PAG data
    pagimg=pagdata.pagimg;
    acid_act_xyz_idx=pagdata.acid_act_xyz_idx;

    %%% polymer distribution matrix
    polym_img=polymdata.polym_img;
    SE_act_xyz_idx=polymdata.SE_act_xyz_idx;

    % pagimg=event.pag.img;
    pag_rcnrad=event.pag.rcnrad;
    pag_grid.x=event.univ.grid.x;
    pag_grid.y=event.univ.grid.y;
    pag_grid.z=event.univ.grid.z;

    pag_density=event.pag.rho; % pag density [/nm3 converted to /cm3];
    Moldensity=1.2/120*6.02*1e23; % molecules/cm3

    while Eold>scatt_Elim
        fprintf(logfile_fid,'............trajcalc3: E = %.4f eV\n',Eold);
        %% 1. Identifying the scattering mechanism 
        
        if strcmp(scattdata.vibr.datasrc,'Khakoo')==1
            vibr_ics=scattdata.vibr.ics;
            vibr_ics(isnan(vibr_ics))=0;
            vibrdata_E=vibr_ics(:,1);
            vibrdata_icsT=sum(vibr_ics(:,2:end),2);
            if Eold>max(vibrdata_E) % no vibrational data here, so go with inelasic
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

%         scattType='HE'; %no vibrational.
        
        %% 2. Calculating the IMFP
        %%%% generate the HE one first, needed whether or not we're LE/HE
        % asking for imfp only
        controlparm.onlyimfp=1;            

        if isfield(scattdata.optical,'inel_dcsdata')
            Elossrand_opt=genrandEloss_OptData(scattdata.optical,Eold,scattdata.optical.inel_dcsdata,controlparm);
        else
            Elossrand_opt=genrandEloss_OptData(scattdata.optical,Eold,controlparm);
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

        imfp_vibr=Inf; % by default, no vibrational scattering
        phi_vibr=0; % temporary, never gets used if the below is commented, so long as imfp_vibr=Inf as above
        theta_vibr=0; % temporary, never gets used if the below is commented, so long as imfp_vibr=Inf as above

        %%%%% comment out the below if/else block if disabling vibrational.
        if imfp_vibr==Inf && imfp_opt==Inf
%             fprintf(logfile_fid,'Warning in trajCalc3: IMFP (Vibr) = IMFP (Optical) = Inf\n');
        end
        rng('shuffle');
        randval=rand;
        imfpratio_vibr_opt=(1/imfp_vibr)/((1/imfp_vibr)+(1/imfp_opt));
        imfp=1/(1/imfp_opt+1/imfp_vibr);            

        %% 3. Scattering Logic
        if randval<imfpratio_vibr_opt % so ifimfp_ratio_vibr_opt=0,never trigger the vibrational path
            %% Vibrational Scattering - obsolete
            Eloss_val=Eloss_vibr;
%             imfp=imfp_vibr; % don't do this
            theta=theta_vibr;
            phi=phi_vibr;
            act='vibr';
        elseif Eold>event.lowEthr
            %% 3.1 Mermin based collision event simulation
            controlparm.onlyimfp=0;
            if isfield(scattdata.optical,'inel_dcsdata')
                Elossrand_opt=genrandEloss_OptData(scattdata.optical,Eold,scattdata.optical.inel_dcsdata,controlparm);
            else
                Elossrand_opt=genrandEloss_OptData(scattdata.optical,Eold,controlparm);
            end
            Eloss_opt=Elossrand_opt.Eloss; % the optical energy loss
            theta_opt=Elossrand_opt.theta; % the optical theta
            phi_opt=Elossrand_opt.phi; % the optical phi

            Eloss_val=Eloss_opt;
            theta = theta_opt;
            phi=phi_opt;

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
        else
            %% 3.2 Low energy random walk (Currently inaccessible)
            act='LowEnergy';
            theta = acos(2*rand-1);
            phi=2*pi*rand;
            
        end
        
        % The angles in the scattering results are relative to the
        % pre-event travelling direction. Thus a coordinate transformation
        % is needed
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Debug: To select a known form of scattering
        %%%================================================================
            scatteringDebug = 'orthogonal';
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
        
        
        %% 4. Post-event electron propagation
        rnew=exprnd(imfp); % exponential distribution        

        znew=zold+rnew*newUnitVec(3);
        xnew=xold+rnew*newUnitVec(1);
        ynew=yold+rnew*newUnitVec(2);

        vec2D_mag=sqrt((xnew-xold)^2+(ynew-yold)^2);
        phival=atan2(ynew-yold,xnew-xold);
        
        %Update the iterators
        theta_old=atan2(vec2D_mag,znew-zold);
        phi_old = phival;

        %%%% calculate the path length
        pathlen=sqrt((xnew-xold)^2+(ynew-yold)^2+(znew-zold)^2);

        %% 5. Energy Deposition [PAG activations? SE-gen?]
  %     [pagidx,npags,polymidx,npolyms]=pag_v_polym([xnew ynew znew],pag_grid,pagimg,polym_img,pag_rcnrad); % move-then-react approach
        [pagidx,npags,polymidx,npolyms]=pag_v_polym([xold yold zold],pag_grid,pagimg,polym_img,pag_rcnrad); % react-then-move approach
        pag_ratio=npags/(npags+npolyms);
        npags_removed=0;

        %% Resolving Low energy electron decisions
        %%% resolve the actual 
        if strcmp(act,'LowEnergy')
            if npags>0
                Eloss_val=pag_Eamin;
                act='LowEnergy-Acid';
            else
                Eloss_val=E_ionize_min;
                act='SE';
            end
        end
        
        if rand<pag_ratio || strcmp(act,'6eVRes') || strcmp(act,'LowEnergy-Acid')% if activating a PAG
         %% 5.1 Acid generation (by volume ratio and 6eV resonance)
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
                num2remove=1; % how many PAGS to remove [could be Eloss/Eact in the future!]
                %%% Determine which pag to remove 
                remove_idx=randi([1 length(pagidx)],num2remove);
                %%% Take the pag out of the image
                pagimg(pagidx(remove_idx))=pagimg(pagidx(remove_idx))-1; % UNCOMMENT TO ENABLE PAG SATURATION EFFECTS
                npags_removed=npags_removed+1;
                nacid=1;
                Ese=Eloss_val-pag_Eamin;    Ese(Ese<0)=0;
                if Ese>0
                    nSE=1;
                    nion=1;
                    SE_act_xyz_idx=[SE_act_xyz_idx pagidx(remove_idx)];
    %                 polym_img(pagidx(remove_idx))=polym_img(pagidx(remove_idx))+1;
                end
                acid_act_xyz_idx=[acid_act_xyz_idx pagidx(remove_idx)];
                act='acid';
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
                    SE_act_xyz_idx=[SE_act_xyz_idx polymidx(add_idx)];
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
                %% 5.2 Secondary electron generation
                    gen_SE=1; % if 1, it will generate secondaries, if 0 then no SE is generated
                    if npolyms>0 % if there was no polymer molecule
                        Ese=(Eloss_val>E_ionize_min).*(Eloss_val-E_ionize_min);
                        nSE=double(Eloss_val>E_ionize_min); % set to 1 if above is set up for SE-gen instead of creating excited states
                        nion=nSE;
                        num2add=nion;
                        if num2add>0 && ~isempty(polymidx)
                            add_idx=randi([1 length(polymidx)],num2add);
%                             polym_img(polymidx(add_idx))=polym_img(polymidx(add_idx))+1;
                            SE_act_xyz_idx=[SE_act_xyz_idx polymidx(add_idx)];
                        end                        
                    else
                        Eloss_val=0; % set to 0 if there was no event
                    end
                case '6eVRes'   %%% !!!! This branch is not reachable
                    act='6eVRes-Polym';
                case 'vibr'     %%% !!!! This branch is not reachable
                    %nacid=1;% temporary, added on 5/27/2017 to test saturation 
                    nacid_unsat=1;
                    if npags>0
        %             Eloss_val=min([6.8 max([Eloss_val pag_Eamin])]); % modify this value, as its a 6.8 eV event instead.
                        num2remove=1; % how many PAGS to remove [could be Eloss/Eact in the future!]
                        remove_idx=randi([1 length(pagidx)],num2remove);
                        pagimg(pagidx(remove_idx))=pagimg(pagidx(remove_idx))-1;  % UNCOMMENT TO ENABLE PAG SATURATION EFFECTS
                        npags_removed=npags_removed+1;
                        nacid=1;
                        nacid_unsat=1;
                        Eloss_val=pag_Eamin;
                        Ese=Eloss_val-pag_Eamin;    Ese(Ese<0)=0;
                        if Ese>0
                            nSE=1;
                            nion=1;
                            SE_act_xyz_idx=[SE_act_xyz_idx pagidx(remove_idx)];
    %                         polym_img(pagidx(remove_idx))=polym_img(pagidx(remove_idx))+1;
                        end
                        acid_act_xyz_idx=[acid_act_xyz_idx pagidx(remove_idx)];
                        act='vibr-acid';
                    else
                        act='vibr-polym';
                    end
            end
        end

        Enew=Eold-Eloss_val;

        fprintf(logfile_fid,'............trajcalc3: [npags,len(pagidx),npags_removed,npolyms,pagratio,nacid,nion,nSE,Eloss,Ese,type] = [%d,%d,%d,%d,%.4f,%d,%d,%d,%.3f,%.3f,%s]\n',npags,length(pagidx),npags_removed,npolyms,pag_ratio,nacid,nion,nSE,Eloss_val,Ese,act);

        xyzglobal.x=[xyzglobal.x xnew];
        xyzglobal.y=[xyzglobal.y ynew];
        xyzglobal.z=[xyzglobal.z znew];
        %%%% post energy-loss variables:
        events{count}.xyz=[xnew ynew znew];
        events{count}.xyzglobal=xyzglobal;
        events{count}.pathlen=pathlen;
        events{count}.Ein=Eold;
        events{count}.Eout=Enew;
        events{count}.Eloss=Eloss_val;
        events{count}.act=act;
        events{count}.imfp=imfp;
        events{count}.rnew=rnew;
        events{count}.theta=theta;
        events{count}.phi=phi;

        events{count}.Ese=Ese;
        events{count}.nacid=nacid;
        events{count}.nacid_unsat=nacid_unsat;
        events{count}.nSE=nSE;
        %%%%% update the pag image
    %     events{count}.pag=event.pag;
    %     events{count}.pag.img=pagimg;
        events{count}.pag.rcnrad=pag_rcnrad;
        events{count}.pag.rho=pag_density;
    %     events{count}.pag.acid_act_xyz_idx=acid_act_xyz_idx;
        events{count}.univ=event.univ;

        %%%%% update the "old" variable values
        Eold=Enew;
        if length(Eold)>1
            dbg=1;
        end
        xold=xnew;
        yold=ynew;
        zold=znew;

        count=count+1;
    end
    pagdata.pagimg=pagimg;
    pagdata.acid_act_xyz_idx=acid_act_xyz_idx;

    polymdata.polym_img=polym_img;
    polymdata.SE_act_xyz_idx=SE_act_xyz_idx;
end
