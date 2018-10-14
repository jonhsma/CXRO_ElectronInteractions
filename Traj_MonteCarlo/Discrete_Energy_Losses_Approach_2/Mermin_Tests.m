%% Mermin function tests: specially, thru-q and thru-w
clear;
clc;
close all;

%%%% The mermin fit parameters:
A_merm=[0.1152 0.4 0.012];
Ei_merm=[6.8 26.5 55];
gamma_merm=[6 13 32];

%%%%% physical parameters:
c=3*1e8; % m/s
ao=53*1e-12;
me=9.11*1e-31; % kg
JeV=1.6*1e-19;
hbar_J=6.626*1e-34/(2*pi); % eV.s
hbar=hbar_J./JeV;
na=1.2/120*6.02*1e23; % molecules/cm3
Ef=10; %eV
vf=sqrt(2*Ef*JeV/me);
qf=sqrt(2*me*Ef*JeV)/hbar_J;

%%%% sweep parameters
Eloss=linspace(0,100,100);
q=linspace(0,100,100).*1e9; % specify as per nm, then scale to /m (SI)

ELFmat=[];

for Eloss_count=1:length(Eloss)
    wval=Eloss(Eloss_count)/hbar;
    for q_count=1:length(q)
        qval=q(q_count);
        ELF=0;
        for i = 1:length(A_merm)
            wi=Ei_merm(i)/hbar;
            gammaval=gamma_merm(i)/hbar;
            epsL=epsMerm(wval,wi,qval,qf,vf,gammaval);
%             epsL=epsMerm(wval,wi,qperp,qf,vf,gammaval);

            ELF_tmp=A_merm(i).*imag(-1./epsL);
            if ~isnan(ELF_tmp)
                ELF=ELF+A_merm(i).*imag(-1./epsL);
            end
        end
        ELFmat(q_count,Eloss_count)=ELF;
    end
end

figure;
imagesc(ELFmat);colorbar;title('ELF');