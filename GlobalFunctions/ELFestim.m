function ELF=ELFestim(E,resist)
JeV=1.6*1e-19;
hbar_J=6.626*1e-34/(2*pi);
hbar=hbar_J/JeV;
m=9.11*1e-31;
Ef=10; %eV
vf=sqrt(2*Ef*JeV/m);
qf=sqrt(2*m*Ef*JeV)/hbar_J;

if strcmp(resist,'fuji')==1
    Ei=[6.8 26.5 55];
    wi=Ei./hbar;
    gamma=[6 13 32];
    Ai=[0.11522 0.4 0.012];
end
if strcmp(resist,'Inpria')==1
    Ei=[5.5 7.2 16.5 21 30 55 85];
    wi=Ei./hbar;
    gamma=[1 1 5 6 20 40 30];
    Ai=[0.005 0.008 0.2 0.17 0.12 0.05 0.008];
end

wi=Ei./hbar;

qo=0;

for i = 1:length(Ei)
    Ep=Ei(i);
%     ELF_D(i,:) = Ai(i).*ELFDrude(Efit,1,Ep,gamma(i),qo,Ef);
    epsL=epsMerm(E./hbar,wi(i),0.1*1e9,qf,vf,gamma(i)/hbar);
    ELF_D(i,:)=Ai(i).*imag(-1./epsL);
end
ELF=sum(ELF_D,1);
