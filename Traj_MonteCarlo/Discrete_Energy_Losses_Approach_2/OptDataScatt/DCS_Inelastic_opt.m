function dsigdOmega=DCS_Inelastic_opt(mermparms,Ei,thetasweep_parms)

% A_merm=[0.1152 0.4 0.012];
% Ei_merm=[6.8 26.5 55];
% gamma_merm=[6 13 32];

%%%% physical parameters
c=3*1e8; % m/s
ao=53*1e-12;
me=9.11*1e-31; % kg
JeV=1.6*1e-19;
hbar_J=6.626*1e-34/(2*pi); % eV.s
hbar=hbar_J./JeV;
na=1.2/120*6.02*1e23; % molecules/cm3
Ef=5; %eV
vf=sqrt(2*Ef*JeV/me);
qf=sqrt(2*me*Ef*JeV)/hbar_J;

A_merm=mermparms.Ai;
Ei_merm=mermparms.Ei;
gamma_merm=mermparms.gamma_i;

diffsig=[];
diffsig_norm=[];

qsq=[];

theta=thetasweep_parms.theta;
delta_theta=thetasweep_parms.delta_theta;

% for count = 1:length(E)
%     Ei=E(count);
    [beta,gamma]=betacalc(Ei);
    v=beta*c;
    ki=sqrt(2*me/((hbar*JeV)^2).*Ei*JeV);
    Eloss=linspace(Ef+1,Ei-1,100);
%     Elossmat(:,count)=Eloss';
    dsigdOmega=[];
    for theta_count=1:length(theta)
        thetaval=theta(theta_count);
        integrand=[];
        for Eloss_count=1:length(Eloss)
            Efinal=Ei-Eloss(Eloss_count);
            kf=sqrt(2*me/((hbar*JeV)^2).*Efinal*JeV);
    %         qsq(count2,:)=kf^2+ki^2-2*ki*kf*cos(theta);
            wval=Eloss(Eloss_count)/hbar; % frequency associated with energy loss value
            qperp=kf*sin(thetaval);
            qpar=ki-kf*cos(thetaval);
            qval=sqrt(qperp^2+qpar^2);
            ELF=0;
            for i = 1:length(A_merm)
                wi=Ei_merm(i)/hbar;
                gammaval=gamma_merm(i)/hbar;
                epsL=epsMerm(wval,wi,qval,qf,vf,gammaval);
    %                 epsL=epslind(wval,wi(l),qvec(k),qf,vf);
                if ~isnan(epsL) % in case epsilon=Delta function, 1/epsilon is 0 
                    ELF=ELF+A_merm(i).*imag(-1./epsL);
                end
            end
%             integrand_term1=2*sin(thetaval)/(pi*ao*me*v^2*na);
            integrand_term1=1/(pi^2*ao*me*v^2*na*1e6); % Turn na into /m3
            integrand_term2=cos(thetaval)*ELF;
            integrand_term3=(sin(thetaval))^2+(Eloss(Eloss_count)*JeV./(hbar_J*kf*v))^2;
            integrand(Eloss_count)=integrand_term1.*integrand_term2./integrand_term3;
        end
        dsigdOmega(theta_count)=trapz(Eloss.*JeV,integrand); % convert Energy to J for correct integration units
    end
% end
dsigdOmega=dsigdOmega*1e4; % Turn into cm2;
sigTot=cumsum(dsigdOmega.*sin(theta).*2*pi).*delta_theta;
