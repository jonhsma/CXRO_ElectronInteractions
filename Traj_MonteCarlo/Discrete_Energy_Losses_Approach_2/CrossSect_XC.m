function outparms=CrossSect_XC(E)

e=1.6*1e-19;
JeV=1.6*1e-19;
eps0=8.85*1e-12;

Eloss=linspace(1,E,100);

term1=-2*pi*e^4./(4*pi*eps0)^2;
term2=1./(E.*JeV);
term3=(1./(Eloss.^2)+1./(E.^2)-1./(Eloss.*(E-Eloss)))./JeV^2;

%%%%% All terms upto now are in SI
dsigdEloss=term1.*term2.*term3.*JeV.*1e4; % convert to cm2/eV

outparms.Eloss=Eloss;
outparms.dsigdEloss=dsigdEloss;