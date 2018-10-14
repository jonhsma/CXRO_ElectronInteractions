function events=trajcalc2(event,scattdata,scatt_Elim,xyzglobal)
pag_Eamin=3;
E_ionize_min=10;
Eold=event.Ese;
xold=event.xyz(1);
yold=event.xyz(2);
zold=event.xyz(3);

Enew=Eold;
count=1;
pagact_prob=0.1;

pagimg=event.pag.img;
pag_rcnrad=event.pag.rcnrad;
pag_grid.x=event.univ.grid.x;
pag_grid.y=event.univ.grid.y;
pag_grid.z=event.univ.grid.z;

pag_density=event.pag.rho; % pag density [/nm3 converted to /cm3];
Moldensity=1.2/120*6.02*1e23; % molecules/cm3

while Eold>scatt_Elim
%     u=rand;
    u=1; % setting to 1 permanenetly disables direct excitation
    if u<=pagact_prob
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
    else
        %%%%%% pick a scattering mechanism
        vibr_ics=scattdata.vibr.ics;
        vibr_ics(isnan(vibr_ics))=0;
        vibrdata_E=vibr_ics(:,1);
        vibrdata_icsT=sum(vibr_ics(:,2:end),2);
        
        if Eold>max(vibrdata_E) % no vibrational data here, so go with inelasic
            scattType='HE'; % scattering is of "high-energy" type
        else
            scattType='LE'; % consider "low-energy" as well
        end
%         scattType='HE'; %no vibrational.
        
        %%%% generate the HE one first, needed whether or not we're LE/HE
        if isfield(scattdata.optical,'inel_dcsdata')
            Elossrand_opt=genrandEloss_OptData(scattdata.optical,Eold,scattdata.optical.inel_dcsdata);
        else
            Elossrand_opt=genrandEloss_OptData(scattdata.optical,Eold);
        end
%         Eloss_opt=Elossrand.Eloss;
        imfp_opt=Elossrand_opt.imfp; % imfp in nm
        imfp_vibr=Inf; % by default, no vibrational scattering
        theta_opt=Elossrand_opt.theta;
        phi_opt=Elossrand_opt.phi;
        
        if strcmp(scattType,'LE')==1
            Elossrand_vibr=genrandEloss_Vibr(scattdata.vibr,Eold);
            Eloss_vibr=Elossrand_vibr.Eloss;
            invimfp=Elossrand_vibr.ics*Moldensity;
            imfp_vibr=1/invimfp*1e7; % imfp in nm
            theta_vibr=Elossrand_vibr.theta;
            phi_vibr=Elossrand_vibr.phi;
        end
        
        %%% assign variables Eloss_val,imfp,theta,phi
        randval=rand;
        imfpratio_vibr_opt=(1/imfp_vibr)/((1/imfp_vibr)+(1/imfp_opt));
        imfp=1/(1/imfp_opt+1/imfp_vibr);
        if randval<imfpratio_vibr_opt
            Eloss_val=Eloss_vibr;
%             imfp=imfp_vibr; % don't do this
            theta=theta_vibr;
            phi=phi_vibr;
            act='vibr';
        else
%             Eloss_val=Eloss_opt;
            Eloss_val=Elossrand_opt.Sp*imfp;
%             imfp=imfp_opt;
            theta=theta_opt;
            phi=phi_opt;
            if Eloss_val==0
                act='none'; % if Sp=0, no energy loss.
            else
                if Eloss_val<E_ionize_min
                    act='6eVRes';
                else
                    if Elossrand_opt.inel_choice==1
                        act='6eVRes';
                    else
                        act='SE';
                    end
                end
            end
        end

    end
    
    %%%%%% calculate co-ordinates [now these parms are output from genrandEloss]
%     imfp=interp1(imfpdata.E,imfpdata.imfp,Eold);
%     theta=pi*rand;
%     phi=-pi+2*pi*rand;
    
    %%%%% the correct one
    rnew=imfp*randn;
    if length(xyzglobal.x)==1 % no previous history present
        znew=zold+rnew*cos(theta);
        xnew=xold+rnew*cos(phi)*sin(theta);
        ynew=yold+rnew*sin(phi)*sin(theta);
    else
        x1=xyzglobal.x(end-1);   x2=xyzglobal.y(end);
        y1=xyzglobal.y(end-1);   y2=xyzglobal.y(end);
        z1=xyzglobal.z(end-1);   z2=xyzglobal.z(end);
        vec2D_mag=sqrt((x2-x1)^2+(y2-y1)^2);
        phival=atan2(y2-y1,x2-x1);
        thetaval=atan2(z2-z1,vec2D_mag);
        
        znew=zold+rnew*cos(theta+thetaval);
        xnew=xold+rnew*cos(phi)*sin(theta+thetaval);
        ynew=yold+rnew*sin(phi)*sin(theta+thetaval);
    end
       
    %%%%% Update the event structure
    event.x=xnew; % update as need this for the PAG deplection effects
    event.y=ynew; % Also, no need to update other variables because this structure is never used, things are stored as variables
    event.z=znew; % But, THESE HERE ARE IMPORTANT
    
    % do this so you have full history, needed to model forward scattering
%     event.x=[event.x xnew];
%     event.y=[event.y ynew];
%     event.z=[event.z znew];
    
%     events{count}.x=xnew;
%     events{count}.y=ynew;
%     events{count}.z=znew;
%     events{count}.theta=theta;
    events{count}.controlparms=event.controlparms;
    
    %%%%%% Decide what to do with the deposted energy [PAG activations?]
    switch act
        case 'SE'
            Ese=Eold-Eloss_val;
            Enew=0;
            nacid=0;
            nSE=1;
        case '6eVRes' % if pag available, pick pag vs. polymer using random number; 
            [yesno,navail,pagidx]=pag_avail([xnew ynew znew],pag_grid,pagimg,pag_rcnrad);
            npolym_molecules=4/3*pi*pag_rcnrad^3*Moldensity*1e-21; % convert Moldensity to /nm3 first
            pagload_ratio=navail/(navail+npolym_molecules);
            
            Ese=0;
            nSE=0;
            
            if yesno==1
                Eloss_val=min([6.8 max([Eloss_val pag_Eamin])]); % modify this value, as its a 6.8 eV event instead.
                remove_idx=randi([1 length(pagidx)],1);
                pagimg(pagidx(remove_idx))=pagimg(pagidx(remove_idx))-1;
                nacid=1;
            else
                Eloss_val=0; % if no pag available, no energy loss.
                nacid=0;
                act='6eVRes-none';
            end
            
            Enew=Eold-Eloss_val;
        otherwise
            nSE=0;
%             Ese=Eold-Eloss_val;
            Ese=0;
            Enew=Eold-Eloss_val;
            nacid=0;
    end
    
    %%%% post energy-loss variables:
    
    events{count}.xyz=[xnew ynew znew];
    events{count}.Ein=Eold;
    events{count}.Eout=Enew;
    events{count}.Eloss=Eloss_val;
    events{count}.act=act;
    events{count}.imfp=imfp;
    events{count}.theta=theta;
    events{count}.phi=phi;
    
    events{count}.Ese=Ese;
    events{count}.nacid=nacid;
    events{count}.nSE=nSE;
    %%%%% update the pag image
%     events{count}.pag=event.pag;
    events{count}.pag.img=pagimg;
    events{count}.pag.rcnrad=pag_rcnrad;
    events{count}.pag.rho=pag_density;
    events{count}.univ=event.univ;
    
    %%%%% update the "old" variable values
    Eold=Enew;
    xold=event.x;
    yold=event.y;
    zold=event.z;
    
    count=count+1;
end
