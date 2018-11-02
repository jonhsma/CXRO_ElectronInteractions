function Elossrand=genrandEloss_OptData(varargin)

    %% Initilization
    % If there is not enough input arguments
    if nargin<2
        error('genrandEloss_OptData: Need at least 2 inputs\n');
    end
    % If there are exactly 2 arguments
    optdata=varargin{1};
    Eo=varargin{2};
    do_dcscalc=1;
    % If there are 3 arguments, DCS doesn't need to be balculated
    if nargin>=3
        dcs_datafile=varargin{3};
        do_dcscalc=0;
    end

    %% Calculating IMFP
    Evec=optdata.E;
    imfpvec=optdata.imfp;
    Spvec=optdata.Sp;

    Efit=linspace(0,max(Evec),1000);
    imfp_fit=interp1(Evec,imfpvec,Efit,'linear','extrap');
    Sp_fit=interp1(Evec,Spvec,Efit,'linear','extrap');
    Sp_fit(Sp_fit<0)=0;

    % imfpval=interp1(Efit,imfp_fit,Eo);
    imfpval=interp1(optdata.E,optdata.imfp,Eo,'linear');
    if isnan(imfpval)
    %     warning('WARNING in genrandEloss_optData: imfpval = NaN; Energy = %.4f eV; Setting imfp to Inf!!\n',Eo);
        imfpval=Inf;
    end
    % Spval=interp1(Efit,Sp_fit,Eo);
    Elossrand.imfp=imfpval;

    if nargin>=4
        controlparms=varargin{4};
        if controlparms.onlyimfp~=0
            return;
        end
    end
    
    %{
    %%%% Pick one of the several inelastic mechanisms [CAREFUL: ASSIGNING Sp
    %%%% Below]
    % inel_choice=choose_inelastic_mech(optdata.merm,Eo);
    % if inel_choice.idx==1
    %     act='6eVRes';
    %     Spval=interp1(optdata.merm{1}.E,optdata.merm{1}.Sp,Eo);
    % else
    %     act='SE';
    %     Spval=interp1(optdata.merm{inel_choice.idx}.E,optdata.merm{inel_choice.idx}.Sp,Eo);
    % end
    % Elossval=imfpval*Spval;
    %}

    if do_dcscalc==1 
        %% Calculate the scattering directly
        % WARNING: If file fails, its likely due to this variable being not set, as I never use this block of code
        %%%% Calculate forward scatt. theta based on conservation rule calculations
        
        theta_sweep.delta_theta=1;
        theta_sweep.theta=[0:theta_sweep.delta_theta:90].*pi/180;

        mermparms.Ai=[0.1152 0.4 0.012];
        mermparms.Ei=[6.8 26.5 55];
        mermparms.gamma_i=[6 13 32];

        thetavec=theta_sweep.theta;
        dsigdOmega=DCS_Inelastic_opt(mermparms,Eo,theta_sweep);
    else
    %     dcs_data=load(dcs_datafile);
        dcs_data=optdata.inel_dcsdata;
        idx=find(dcs_data.E(1,:)==Eo);
        dsigdOmega_fit=[];
        if ~isempty(idx)
    %         thetavec=dcs_data.theta; % just a single vector in the most recent one
    %         dsigdOmega=dcs_data.dsigdOmega(idx,:); % row corresponds to incident energy

            thetavec=dcs_data.thetamat(idx,:); % just a single vector in the most recent one
            dsigdOmega=dcs_data.dsigdOmega(idx,:); % row corresponds to incident energy

    %         thetavec2=linspace(1e-3,pi/2,100);
    %         dsigdOmega=interp1(thetavec,dsigdOmega,thetavec2);
    %         thetavec=thetavec2;

            Elossvec=dcs_data.Elossmat(idx,:);
            dsigdE=dcs_data.dsigdE(idx,:);
        else
            dsigdOmega1=dcs_data.dsigdOmega;
            dsigdE1=dcs_data.dsigdE;
            idx1=find(dcs_data.E(1,:)<Eo);

            if isempty(idx1)
    %             thetavec=dcs_data.theta;
    %             dsigdOmega=dcs_data.dsigdOmega(1,:);

                thetavec=dcs_data.thetamat(1,:);
                dsigdOmega=dcs_data.dsigdOmega(1,:);

                Elossvec=dcs_data.Elossmat(1,:);
                dsigdE=dcs_data.dsigdE(1,:);
            else
                idx1=idx1(end);
                idx2=find(dcs_data.E(1,:)>Eo);
                idx2=idx2(1); % energy always >92, so this shouldn't be a problem by default
                dsigdOmega=[];
                dsigdE=[];
    %             thetavec=dcs_data.theta;
                thetavec=[];
                Elossvec=[];

                %%%% fit the dsigdOmega
                theta1=dcs_data.thetamat(idx1,:);
                dcs1=dcs_data.dsigdOmega(idx1,:);
                theta2=dcs_data.thetamat(idx2,:);
                dcs2=dcs_data.dsigdOmega(idx2,:);

                thetavec_min=1e-3;
                thetavec=linspace(thetavec_min,pi/2,100);
                dcs1B=interp1(theta1,dcs1,thetavec,'linear','extrap');
                dcs2B=interp1(theta2,dcs2,thetavec,'linear','extrap');

                for i = 1:size(dcs1B,2)
                    p=polyfit(dcs_data.E(1,idx1:idx2),[dcs1B(i) dcs2B(i)],1);
                    dsigdOmega(1,i)=polyval(p,Eo);
                end

                dbg=1;
                %%%% fit the dsigdE
                Eloss1=dcs_data.Elossmat(idx1,:);
                dcs1=dcs_data.dsigdE(idx1,:);
                dcs1=dcs1(Eloss1~=0);
                Eloss1=Eloss1(Eloss1~=0);

                Eloss2=dcs_data.Elossmat(idx2,:);
                dcs2=dcs_data.dsigdE(idx2,:);
                dcs2=dcs2(Eloss2~=0);
                Eloss2=Eloss2(Eloss2~=0);

                delta_E=0.1;
                Elossmin=min([min(Eloss1) min(Eloss2)]);
                Elossmax=Eo;
                Elossvec=Elossmin:delta_E:Elossmax;
                if length(Elossvec)<=3
                    Elossvec=linspace(Elossmin,Elossmax,3);
                end

    %             Elossvec=linspace(min([min(Eloss1) min(Eloss2)]),Eo,100);
    %             try
                    dcs1B=interp1(Eloss1,dcs1,Elossvec,'linear','extrap');
    %             catch errorstr
    %                 dbg=1;
    %             end

    %             try
                    dcs2B=interp1(Eloss2,dcs2,Elossvec,'linear','extrap');
    %             catch errorstr
    %                 dbg=1;
    %             end

                dcs1B(Elossvec>max(Eloss1))=0;
                dcs2B(Elossvec>max(Eloss2))=0;

                for i = 1:size(dcs1B,2)
                    %%%% fit the dsigdE
    %                 p=polyfit(dcs_data.E(1,idx1:idx2)',dcs_data.dsigdE(idx1:idx2,i),1);
    %                 dsigdE(1,i)=polyval(p,Eo);

                    p=polyfit(dcs_data.E(1,idx1:idx2),[dcs1B(i) dcs2B(i)],1);
                    dsigdE(1,i)=polyval(p,Eo);
                end
                dbg=1;
            end
        end
    end

    dsigdE(dsigdE<0)=0;
    % if max(Elossvec)<Eo
    %     Elossvec=[Elossvec Eo];
    %     dsigdE=[dsigdE 0];
    % end

    if any (dsigdOmega<0)
        dbg=1;
        dsigdOmega(dsigdOmega<0)=0; % happens if interpolating out to 90 degrees, becomes a small negative number. JUST BE AWARE FOR THE FUTURE
    end
    %%%% sample theta and phi randomly
    % dcs_cdf=cumsum(dsigdOmega);
    
    %% Calculated the culmulative distribution of the differential cross-section
    dcs_cdf     =   [];
    dcs_cdf(1)  =   0;
    for i = 2:length(dsigdOmega)
        dcs_cdf(i)=trapz(thetavec(1:i),dsigdOmega(1:i));
        %dcs_cdf(i)=trapz(thetavec(1:i),sin(thetavec(1:i)).*dsigdOmega(1:i));
    end
    dcs_cdf=dcs_cdf-dcs_cdf(1);
    dcs_cdf=dcs_cdf./dcs_cdf(end);

    theta_rand=randgen(thetavec,dcs_cdf,1);
    phi_rand=2*pi*rand;

    %%%% sample Energy loss
    dcs_cdf=cumsum(dsigdE);
    dcs_cdf=[];
    for i = 2:length(dsigdE)
        dcs_cdf(i)=trapz(Elossvec(1:i),dsigdE(1:i));
    end
    dcs_cdf=dcs_cdf-dcs_cdf(1);
    dcs_cdf=dcs_cdf./dcs_cdf(end);

    Eloss_rand=randgen(Elossvec,dcs_cdf,1);
    % Eloss_rand=trapz(Elossvec,Elossvec.*dsigdE)/trapz(Elossvec,dsigdE); % Disable randomness for debugging purposes.
    % maxidx=find(dsigdE==max(dsigdE(:)));
    % Eloss_rand=Elossvec(maxidx);

    % if max(Elossvec)<=100
    %     dbg=1;
    % end

    %%%% sample cross-section
    % if exist('idx')
    %     E_TCS=dcs_data.icsdata{idx}.E;
    %     theta_TCS=dcs_data.icsdata{idx}.theta;
    %     TCS_2D=dcs_data.icsdata{idx}.CS_2D;
    %     dbg=1;
    % else
    % end

    Elossrand.Eloss=Eloss_rand;
    % Elossrand.inel_choice=inel_choice.idx;
    % Elossrand.Sp=Spval;
    % if Eo>40
    %     Elossrand.imfp=1.2*imfpval;
    % else
    % end

    rng('shuffle');
    if rand<=0.5
        theta_rand=-1*theta_rand; % equal probaiblity for [-pi/2,0], and [0,pi/2]
    end

    % theta_rand=-pi+pi*rand; % temporary just to see.

    Elossrand.theta=theta_rand;
    % Elossrand.theta=pi*rand; % uncomment only when you want to randomize the angle 
    Elossrand.phi=phi_rand;
end
