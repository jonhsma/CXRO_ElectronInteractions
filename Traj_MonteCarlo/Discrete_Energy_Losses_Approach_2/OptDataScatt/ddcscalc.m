function outparms=ddcscalc(varargin)

if nargin<1
    error('ddcs: Need at least 1 argument, the incident electron energy\n');
end

E=varargin{1};
%%%% The mermin fit parameters: OLD
A_merm=[0.1152 0.4 0.012];
Ei_merm=[6.8 26.5 55];
gamma_merm=[6 13 32];

%%%% The mermin fit parameters: Newest, lsqfit 
A_merm=[0.051 0.295 0.224 0.001];
Ei_merm=[6.22 23.22 32.58 295.7];
gamma_merm=[3.14 10.69 26.49 24];

% A_merm=[0.4 0.012];
% Ei_merm=[26.5 55];
% gamma_merm=[13 32];
    
if nargin==2
    A_merm=mermparms.Ai;
    Ei_merm=mermparms.Ei;
    gamma_merm=mermparms.gamma;
end

%%%%% physical parameters:
c=3*1e8; % m/s
ao=53*1e-12;
me=9.11*1e-31; % kg
JeV=1.6*1e-19;
hbar_J=6.626*1e-34/(2*pi); % eV.s
hbar=hbar_J./JeV;
na=1.2/120*6.02*1e23; % molecules/cm3

Ef=15.5; %eV
Eg=1e-3; % so that the logarithmic Eloss matrix doesn't do werid stuff

vf=sqrt(2*Ef*JeV/me); % m/s
qf=sqrt(2*me*Ef*JeV)/hbar_J; %1/m

%%%% sweep values
theta_min=1e-9;
theta_max=90;
delta_theta=0.5;
theta_d=theta_min:delta_theta:theta_max; % degrees

Esweep_logspace_thr=1000;
if E>Esweep_logspace_thr % large angles, use logarithmic scale for precision
    theta_min=1e-3; % can't be zero for logrithmic scale;
    theta_d=logspace(log10(theta_min),log10(theta_max),length(theta_d));
end
theta=theta_d.*pi/180; % radians

for count = 1:length(E)
    Ei=E(count); % eV
    [beta,gamma]=betacalc(Ei);
    v=beta*c; % m/s
    ki=sqrt(2*me/((hbar*JeV)^2).*Ei*JeV); % 1/m
    
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
    
    dsigdOmega=[];
    dsigdOmega_Sp=[];
    Sp_integrand_theta=[];
    for theta_count=1:length(theta)
        thetaval=theta(theta_count);
                
        integrand=[];
        Sp_integrand_Eloss=[];
        for Eloss_count=1:length(Eloss)
            Efinal=Ei-Eloss(Eloss_count);
            kf=sqrt(2*me/((hbar*JeV)^2).*Efinal*JeV); % 1/m
            wval=Eloss(Eloss_count)/hbar; % hbar is in eV.s unit
            
            %%% set up the sine and cosine values for the theta ranges
            if thetaval<=pi/2
                sineval=sin(thetaval);
                cosval=cos(thetaval);
                qperp=kf*sineval;
                qpar=ki-kf*cosval;
            else
                sineval=sin(pi-thetaval);
                cosval=cos(pi-thetaval);
                qperp=kf*sineval;
                qpar=kf*cosval+ki;
                dbg=1;
            end
            
            qval=sqrt(qperp^2+qpar^2);
            
            ELF=0;
            for i = 1:length(A_merm)
                wi=Ei_merm(i)/hbar; % hbar unit is eV.s
                gammaval=gamma_merm(i)/hbar;
                epsL=epsMerm(wval,wi,qval,qf,vf,gammaval);
%                 epsL=epsMerm(wval,wi,qperp,qf,vf,gammaval);
                
                ELF_tmp=A_merm(i).*imag(-1./epsL);
                
                if ELF_tmp<0
                    ELF_tmp=0;
                end
                
                if ~isnan(ELF_tmp)
                    ELF=ELF+ELF_tmp;
                end
            end
            
            integrand_term1=1/(pi^2*ao*me*v^2*na*1e6); % Turn na into /m3
            
%             if thetaval<=pi/2
%                 integrand_term2=cos(thetaval)*ELF;
%             else
%                 integrand_term2=cos(pi-thetaval)*ELF;
%             end
            integrand_term2=(cosval)*ELF; % Looks NON-PHYSICAL; if theat>pi/2, use symmetry of cosine by using abs(cos(theta))
%             integrand_term2=ELF;
            
            integrand_term3=(sineval)^2+(Eloss(Eloss_count)*JeV./(hbar_J*kf*v))^2;
%             integrand_term3=(sin(thetaval))^2+(Eloss(Eloss_count)*JeV./(kf*v))^2;
%             integrand_term3=(sin(thetaval))^2+(Eloss(Eloss_count)*JeV./(me*v^2))^2;
%             integrand_term3=((thetaval))^2+(Eloss(Eloss_count)*JeV./(hbar_J*kf*v))^2;

            Sp_integrand1=2/(pi*ao*me*v^2);
            Sp_integrand2=cos(thetaval)*ELF*Eloss(Eloss_count)*JeV;
            Sp_integrand3=sin(thetaval)./((sin(thetaval))^2+(Eloss(Eloss_count)*JeV/(hbar_J*kf*v))^2);
            Sp_integrand_Eloss(1,Eloss_count)=Sp_integrand1.*Sp_integrand2.*Sp_integrand3;
            
%             integrand(Eloss_count)=integrand_term1.*integrand_term2./integrand_term3;
%             ddcs(theta_count,Eloss_count,count)=integrand_term1.*integrand_term2./integrand_term3;
%             ddcs(theta_count,Eloss_count,count)=1./(pi^2*ao*me*v^2*na*1e6) .* ELF.*ki.*kf./(ki.^2+kf.^2-2.*ki.*kf.*cosval); % Alternative, not doing omega/velocity from Ritchie
            ddcs(theta_count,Eloss_count,count)=1./(2*pi^2*ao*E) .* ELF.*ki.*kf./(ki.^2+kf.^2-2.*ki.*kf.*cosval); % The correct one, directly using Pine's equation
            integrand(Eloss_count)=ddcs(theta_count,Eloss_count,count);
            
%             ddcs(theta_count,Eloss_count,count)=integrand_term2./integrand_term3;
            ELFval(theta_count,Eloss_count,count)=ELF;
            qvec(theta_count,Eloss_count,count)=qval;
        end
        Sp_integrand_theta(1,theta_count)=trapz(Eloss.*JeV,Sp_integrand_Eloss);
        
        dsigdOmega(theta_count)=trapz(Eloss,integrand); % integrand is /m/eV
        dsigdOmega_Sp(theta_count)=trapz(Eloss.*JeV,Eloss.*integrand); % integrand in /J, so Eloss can be in eV here
    end
end

%%%% calculate Sp through integration of DDCS just to double-check
% Sp_v2=[];
% for i = 1:size(ddcs,2)
%     tmpvec=ddcs(:,i);
%     Sp_v2(i)=trapz(theta',tmpvec.*sin(theta').*2*pi);
% end

% CS=[];
% CS2=[];
% for theta_count=1:length(theta)
%     ddcs_vec=ddcs(theta_count,:);
% %     ddcs_sum=cumsum(ddcs_vec);
% %     CS(theta_count,:)=ddcs_sum(2:end).*(Eloss(2:end)-Eloss(1:end-1)).*JeV;
%     [defint CS(theta_count,:)]=integrate(Eloss.*JeV,ddcs_vec);
% end

% for Eloss_count=1:size(CS,2)
%     ddcs_vec=CS(:,Eloss_count);
% %     ddcs_sum=cumsum(ddcs_vec.*sin(theta').*2*pi);
% %     CS2(:,Eloss_count)=ddcs_sum(2:end).*(theta(2:end)'-theta(1:end-1)');
%     [defint CS2(:,Eloss_count)]=integrate(theta',ddcs_vec.*sin(theta').*2*pi);
% end
% CS2=CS2.*1e4; % 2-D Integrated Cross-Section

% Sp_v2=trapz(theta,Sp_integrand_theta); % J/m
% Sp_v2=Sp_v2*1e-9/(1.6*1e-19); % eV/nm;

% dsigdOmega=dsigdOmega*1e4; % In new version, its DIIMFP, so no need to do cm2 conversion
dsigdOmega_Sp=dsigdOmega_Sp*1e4.*na; % Turn into eV/cm
dsigdOmega_Sp=dsigdOmega_Sp.*1e-7; % convert to eV/nm

ics=cumsum(dsigdOmega.*sin(theta).*2*pi).*delta_theta*pi/180;
ics_Sp=cumsum(dsigdOmega_Sp.*sin(theta).*2*pi).*delta_theta*pi/180; % eV/cm as multiplied by na above already

sigTot=trapz(theta,dsigdOmega.*sin(theta).*2*pi); % sigTot is now inv-IMFP
Sp=trapz(theta,dsigdOmega_Sp.*sin(theta).*2*pi);

% imfp=1./(sigTot.*na).*1e7; % convert to nm
imfp = 1/sigTot.*1e9; % in new version, sigTot is 1/m, so imfp is now nm

% outparms.CrossSect.E=Eloss(2:end);
% outparms.CrossSect.theta=theta(2:end);
% outparms.CrossSect.CS_2D=CS2;

outparms.dsigdOmega=dsigdOmega;
outparms.dsigdOmega_Sp=dsigdOmega_Sp;

outparms.ddcs=ddcs;
outparms.ELF=ELFval;
outparms.qvec=qvec;
outparms.ics=ics;
outparms.imfp=imfp;
outparms.sigTot=sigTot;
outparms.Sp=Sp;
% outparms.Sp2=Sp_v2;
outparms.Eloss=Eloss;
outparms.theta=theta;