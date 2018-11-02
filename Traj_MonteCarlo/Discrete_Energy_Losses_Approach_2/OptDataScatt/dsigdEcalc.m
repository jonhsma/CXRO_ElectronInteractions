function outparms=dsigdEcalc(varargin)

if nargin<1
    error('ddcs: Need at least 1 argument, the incident electron energy\n');
end

E=varargin{1};

%%%% The mermin fit parameters: OLD
%{
A_merm=[0.1152 0.4 0.012];
Ei_merm=[6.8 26.5 55];
gamma_merm=[6 13 32];
%}

%%%% The mermin fit parameters: Newest, lsqfit 
A_merm=[0.051 0.295 0.224 0.001];
Ei_merm=[6.22 23.22 32.58 295.7];
gamma_merm=[3.14 10.69 26.49 24];

% A_merm=[0.4 0.012];
% Ei_merm=[26.5 55];
% gamma_merm=[13 32];
    
if nargin==2
    mermparams     =   varargin{2};
    A_merm          =   mermparams.Ai;
    Ei_merm         =   mermparams.Ei;
    gamma_merm      =   mermparams.gamma_i;
end

%%%%% physical parameters:
c=3*1e8; % m/s
ao=53*1e-12;
me=9.11*1e-31; % kg
JeV=1.6*1e-19;
hbar_J=6.626*1e-34/(2*pi); % J.s
hbar=hbar_J./JeV; % eV.s
na=1.2/120*6.02*1e23; % molecules/cm3

Ef=15.5; %eV
Eg=0.001;

vf=sqrt(2*Ef*JeV/me);
qf=sqrt(2*me*Ef*JeV)/hbar_J;

diffsig=[];
diffsig_norm=[];

qsq=[];

%%%% sweep values
theta_min=0;
theta_max=90;
delta_theta=0.5;
numthetavals=(theta_max-theta_min)/delta_theta;
theta_d=theta_min:delta_theta:theta_max;

Esweep_logspace_thr=200;
if E>Esweep_logspace_thr % large angles, use logarithmic scale for precision
    theta_min=1e-3; % can't be zero for logrithmic scale;
    theta_d=logspace(log10(theta_min),log10(theta_max),length(theta_d));
end
theta=theta_d.*pi/180;

imfp=[];
for count = 1:length(E)
%     fprintf('Energy %d of %d\n',count,length(E));
    Ei=E(count);
    [beta,gamma]=betacalc(Ei);
    v=beta*c;
    ki=sqrt(2*me/((hbar*JeV)^2).*Ei*JeV);
    
%     Eloss_min=1e-3;
    Eloss_min=Eg;
    Eloss_max=Ei-Ef;
%     Eloss_max=Ei;
    
    if E>Esweep_logspace_thr
        Eloss=logspace(log10(Eloss_min),log10(Eloss_max),1000);
    else
        Eloss=linspace(Eloss_min,Eloss_max,500);
%         Eloss=Eloss_min:Eloss_delta:Ei;
    end
%     Elossmat(:,count)=Eloss';
    dsigdOmega=[];
    dsigdOmega_Sp=[];
    Sp_integrand_theta=[];
        
    for i = 1:length(Eloss)
%         fprintf('...Eloss %d of %d\n',i,length(Eloss));
        if i==length(Eloss)
            dbg=1;
        end
        qmin = sqrt(2*me*Ei*JeV)/hbar_J - sqrt(2*me*(Ei-Eloss(i))*JeV)/hbar_J;
        qmax = sqrt(2*me*Ei*JeV)/hbar_J + sqrt(2*me*(Ei-Eloss(i))*JeV)/hbar_J;
        qvec=logspace(log10(qmin),log10(qmax),200);
%         qvec=linspace(qmin,qmax,500);
        wval=Eloss(i)/hbar; % frequency associated with energy loss value
        
        ELF=zeros(size(qvec));
        for j = 1:length(qvec)
            
            qval=qvec(j);
            for k = 1:length(A_merm)
                wi=Ei_merm(k)/hbar;
                gammaval=gamma_merm(k)/hbar;
                epsL=epsMerm(wval,wi,qval,qf,vf,gammaval);
%                 epsL=epsMerm(wval,wi,qperp,qf,vf,gammaval);
                
                ELF_tmp=A_merm(k).*imag(-1./epsL);
                
                if ELF_tmp<0 | ELF_tmp==NaN
                    ELF_tmp=0;
                end
                
%                 if ~isnan(ELF_tmp)
                    ELF(j)=ELF(j)+ELF_tmp;
%                 end
            end
        end
        integrand=1./(pi*ao.*Ei).*ELF.*1./qvec;
        dsigdE(count,i)=trapz(qvec,integrand);
        eLossMat(count,i)=Eloss(i);
    end
    iimfp=trapz(Eloss,dsigdE(count,:));
    imfp=1./iimfp.*1e9; % invert, then convert to nm
end

outparms.eLossMat=eLossMat;
outparms.dsigdE=dsigdE;
outparms.imfp=imfp;
% imfp=1./(sigTot.*na).*1e7; % convert to nm

% outparms.CrossSect.E=Eloss(2:end);
% outparms.CrossSect.theta=theta(2:end);
% outparms.CrossSect.CS_2D=CS2;

% outparms.dsigdOmega=dsigdOmega;
% outparms.dsigdOmega_Sp=dsigdOmega_Sp;

% outparms.ddcs=ddcs;
% outparms.ELF=ELFval;
% outparms.qvec=qvec;
% outparms.ics=ics;
% outparms.imfp=imfp;
% outparms.sigTot=sigTot;
% outparms.Sp=Sp;
% % outparms.Sp2=Sp_v2;
% outparms.Eloss=Eloss;
% outparms.theta=theta;