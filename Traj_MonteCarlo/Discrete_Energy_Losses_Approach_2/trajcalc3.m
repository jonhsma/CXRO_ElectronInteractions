%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What it does: start a trajectory with the secondary from the input event
%%%     scatter -> propagate -> scatter -> ..... and so on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [events,pagdata,polymdata]=trajcalc3(event,scattdata,scatt_Elim,xyzglobal,pagdata,polymdata,logfile_fid)
global illustration;

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
    %     fprintf('...Eold = %.2f eV\n',Eold);
    %     u=rand;
        fprintf(logfile_fid,'............trajcalc3: E = %.4f eV\n',Eold);
        u=1; % setting to 1 permanenetly disables direct excitation
        if u<=pagact_prob % always ignored now, by setting u = 1
            %{
            Ese=0;
            nSE=0;
            if Eold>pag_Eamin
                act='pag';

                Eloss_val=pag_Eamin;
                if event.controlparms.model_pag_devel~=0
                    [event,nacid]=PAG_depl(event,Eloss_val,pag_Eamin);
                else
                    nacid=floor(Eloss_val/pag_Eamin);
                end

    %             nacid=1;

                if nacid>0
                    Enew=Eold-Eloss_val;
                end
            else
                act='none';
                nacid=0;
                Eloss_val=0;
                Enew=Eold;
            end
            %}
        else
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
    %{
    %         Eloss_opt=Elossrand_opt.Eloss; % above is only for getting imfp, so Eloss=0 for above lines in genrandEloss_OptData
    % % %         Evec=scattdata.optical.inel_dcsdata.E;
    % % %         idx=find(Evec==Eloss_opt);
    % % %         if ~isempty(idx)
    % % %             E_TCS=scattdata.optical.CrossSect.E;
    % % %             theta_TCCS=scattdata.optical.CrossSect.theta;
    % % %             TCS_2D=scattdata.optical.CrossSect.CS_2D;
    % % %             dbg=1;
    % % %         else
    % % %             idx1=find(Evec<Eloss_opt);
    % % %         end
    %}

            imfp_vibr=Inf; % by default, no vibrational scattering
            phi_vibr=0; % temporary, never gets used if the below is commented, so long as imfp_vibr=Inf as above
            theta_vibr=0; % temporary, never gets used if the below is commented, so long as imfp_vibr=Inf as above

            %%%%% comment out the below if/else block if disabling vibrational.
    %{
    % % %         if strcmp(scattType,'LE')==1 % Khakoo data from above if-else
    % % % %         if strcmp(scattdata.vibr.datasrc,'Khakoo')==1
    % % %             Elossrand_vibr=genrandEloss_Vibr(scattdata.vibr,Eold);
    % % %             Eloss_vibr=Elossrand_vibr.Eloss;
    % % %             invimfp=Elossrand_vibr.ics*Moldensity;
    % % %             imfp_vibr=1/invimfp*1e7; % imfp in nm
    % % %             theta_vibr=Elossrand_vibr.theta;
    % % %             phi_vibr=Elossrand_vibr.phi;
    % % %         else % Frohlich data from above if-else
    % % %             if strcmp(scattdata.vibr.datasrc,'Frohlich')==1
    % % %                 if Eold<=Inf % just temporary to speed things up [Inf if you want to execute whats in the if-block]
    % % %                     eps0=scattdata.vibr.eps0;
    % % %                     epsinf=scattdata.vibr.epsinf;
    % % %                     E_optphonon=scattdata.vibr.hbarw;
    % % %                     imfp_vibr1=[];
    % % %                     for Eph_count=1:length(E_optphonon)
    % % %     %                     [imfp_vibr,imfp_vibr_creation,imfp_vibr_annihilation,theta_vibr]=scattdata.vibr.imfp_func(Eold,eps0,epsinf,E_optphonon);
    % % %                         [imfp_vibr1(Eph_count),imfp_vibr_creation,imfp_vibr_annihilation,theta_vibr]=scattdata.vibr.imfp_func(Eold,eps0,epsinf,E_optphonon(Eph_count));
    % % %     %                     [imfp_vibr2,imfp_vibr_creation,imfp_vibr_annihilation,theta_vibr]=scattdata.vibr.imfp_func(Eold,eps0,epsinf,E_optphonon(Eph_count));
    % % %                     end
    % % %                     imfp_vibr=1/(sum(1./imfp_vibr1));
    % % %                     dbg=1;
    % % %     %                 imfp_vibr=1/(1/imfp_vibr1 + 1/imfp_vibr2);
    % % %     %                 if rand<(1/imfp_vibr1)/(1/imfp_vibr1 + 1/imfp_vibr2)
    % % %     %                     Eloss_val=E_optphonon(1);
    % % %     %                 end
    % % %                     Eloss_vibr=E_optphonon(1);
    % % %         %             Eloss_vec=scattdata.vibr.E;
    % % %         %             imfp_vec=scattdata.vibr.imfp;
    % % %         %             imfp_vibr=interp1(Eloss_vec,imfp_vec,Eold,'linear','extrap');
    % % %                     theta_vibr=rand*2*pi;
    % % %                     rng('shuffle');
    % % %                     phi_vibr=-pi+2*pi*rand;
    % % %                 else
    % % %                     imfp_vibr=Inf;
    % % %                     phi_vibr=0;
    % % %                     theta_vibr=0;
    % % %                 end
    % % %             end
    % % %         end
            %}
            %%%% the above commented if disabling vibrations
            %%% assign variables Eloss_val,imfp,theta,phi
            if imfp_vibr==Inf & imfp_opt==Inf
    %             fprintf(logfile_fid,'Warning in trajCalc3: IMFP (Vibr) = IMFP (Optical) = Inf\n');
            end
            rng('shuffle');
            randval=rand;
            imfpratio_vibr_opt=(1/imfp_vibr)/((1/imfp_vibr)+(1/imfp_opt));
            imfp=1/(1/imfp_opt+1/imfp_vibr);            
   
            %% 3. Calculating the energy loss and determine the type of action
            if randval<imfpratio_vibr_opt % so ifimfp_ratio_vibr_opt=0,never trigger the vibrational path
                Eloss_val=Eloss_vibr;
    %             imfp=imfp_vibr; % don't do this
                theta=theta_vibr;
                phi=phi_vibr;
                act='vibr';
            else % if not vibrational, extract the ionization IMFP
                if Eold>event.lowEthr
                    controlparm.onlyimfp=0;
                    if isfield(scattdata.optical,'inel_dcsdata')
                        Elossrand_opt=genrandEloss_OptData(scattdata.optical,Eold,scattdata.optical.inel_dcsdata,controlparm);
                    else
                        Elossrand_opt=genrandEloss_OptData(scattdata.optical,Eold,controlparm);
                    end
                    Eloss_opt=Elossrand_opt.Eloss; % the optical energy loss
                    theta_opt=Elossrand_opt.theta; % the optical theta
                    phi_opt=Elossrand_opt.phi; % the optical phi

        %             Eloss_val=Elossrand_opt.Sp*imfp;
                    Eloss_val=Eloss_opt;
        %             imfp=imfp_opt;
                    theta=theta_opt;
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
                    act='LowEnergy';
                    theta=-pi+2*pi*rand;
                    %theta = acos(2*rand-1);
                    phi=2*pi*rand;
                end
            end
        end
        
        %% 4. Post-event electron propagation
        rnew=exprnd(imfp); % exponential distribution
%{
    %     theta_new=theta_old+theta;
    %     if theta_old>pi/2
    % %         theta_new=pi-theta;
    %         cosval=cos(pi-theta);
    %     else
    % %         theta_new=theta;
    %         cosval=cos(theta);
    %     end
%}
        cosval=cos(theta_old+theta);
        sineval=sqrt(1-cosval^2);

        znew=zold+rnew*cosval;
        xnew=xold+rnew*cos(phi)*sineval;
        ynew=yold+rnew*sin(phi)*sineval;

        vec2D_mag=sqrt((xnew-xold)^2+(ynew-yold)^2);
        phival=atan2(ynew-yold,xnew-xold);
        theta_old=atan2(znew-zold,vec2D_mag);

        %%%% calculate the path length
        pathlen=sqrt((xnew-xold)^2+(ynew-yold)^2+(znew-zold)^2);

        %% 5. Energy Deposition [PAG activations? SE-gen?]
  %     [pagidx,npags,polymidx,npolyms]=pag_v_polym([xnew ynew znew],pag_grid,pagimg,polym_img,pag_rcnrad); % move-then-react approach
        [pagidx,npags,polymidx,npolyms]=pag_v_polym([xold yold zold],pag_grid,pagimg,polym_img,pag_rcnrad); % react-then-move approach
        pag_ratio=npags/(npags+npolyms);
        npags_removed=0;

        %% 5.1 Low energy (incident electron) event
        %%%% uncomment below if simulating lower energies
        if strcmp(act,'LowEnergy')
            if npags>0
                Eloss_val=pag_Eamin;
                act='LowEnergy-Acid';
            else
                Eloss_val=E_ionize_min;
                act='SE';
            end
        end
        
        if rand<pag_ratio || strcmp(act,'6eVRes') % strcmp(act,'LowEnergy-Acid')% if activating a PAG
         %% acid generation (by volume ratio and 6eV resonance)
            Ese=0;
            nSE=0;
            nion=0;
            nacid=0;
            nacid_unsat=1;
    %         nacid=1; % temporary, added on 5/27/2017 to test saturation 
            if npags>0
    %             Eloss_val=min([6.8 max([Eloss_val pag_Eamin])]); % modify this value, as its a 6.8 eV event instead.
                num2remove=1; % how many PAGS to remove [could be Eloss/Eact in the future!]
                remove_idx=randi([1 length(pagidx)],num2remove);
                pagimg(pagidx(remove_idx))=pagimg(pagidx(remove_idx))-1; % UNCOMMENT TO ENABLE PAG SATURATION EFFECTS
                npags_removed=npags_removed+1;
                nacid=1;
    %             Eloss_val=pag_Eamin;
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
                    nSE=1;
                    nion=1;
                    num2add=nion;
                    add_idx=randi([1 length(polymidx)],num2add);
                    SE_act_xyz_idx=[SE_act_xyz_idx polymidx(add_idx)];
    %                 polym_img(polymidx(add_idx))=polym_img(polymidx(add_idx))+1;
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
                %% Secondary electron generation
                    gen_SE=1; % if 1, it will generate secondaries, if 0 then no SE is generated
                    if npolyms>0 % if there was no polymer molecule
                        if gen_SE==1
                            Ese=(Eloss_val>E_ionize_min).*(Eloss_val-E_ionize_min);
                            nSE=double(Eloss_val>E_ionize_min); % set to 1 if above is set up for SE-gen instead of creating excited states
                            nion=nSE;
                            num2add=nion;
                            if num2add>0 && ~isempty(polymidx)
                                add_idx=randi([1 length(polymidx)],num2add);
    %                             polym_img(polymidx(add_idx))=polym_img(polymidx(add_idx))+1;
                                SE_act_xyz_idx=[SE_act_xyz_idx polymidx(add_idx)];
                            end
                        end
                    else
                        Eloss_val=0; % set to 0 if there was no event
                    end
                case '6eVRes'
                    %%% !!!! This branch is not reachable
                    act='6eVRes-Polym';
                case 'vibr'
                    %%% !!!! This branch is not reachable
  %                 nacid=1;% temporary, added on 5/27/2017 to test saturation 
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
