function ELF = ELFDrude(E,oscmax,Ep,gamma,q,Ef)
%%% ELF = ELFDrude(E,Ep,gamma,q,Ef)
%%% E = Energy (eV)
%%% oscmax = Max value of the oscillator at Ep
%%% Ep = Plasmon Energy (eV)
%%% gamma = FWHM (eV)
%%% q = Momentum Transfered (/m)
%%% Ef = Fermi Energy (eV)
m=9.11*1e-31;
JeV=1.6*1e-19;
hbar_J=6.626*1e-34/(2*pi);

if Ep~=0
    alpha=3/5*Ef/Ep;
    Epq=Ep+alpha*hbar_J^2/m*q^2/JeV;
    ELF=E.*gamma.*Epq^2./((E.^2-Epq^2).^2+(E.*gamma).^2);
    ELF(ELF<0)=eps;
else
    sig=gamma/2.35;
    gamma=gamma/2;
    ELF=1/(sig*sqrt(2*pi)).*exp(-(E-Ep).^2./(2*sig^2));
%     ELF=1/pi.*gamma/2./((E-Ep).^2+(gamma/2)^2);
end
% ELF=ELF./max(ELF).*oscmax;