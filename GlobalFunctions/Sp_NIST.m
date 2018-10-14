function Sp=Sp_NIST(E,I,Z,Molmass,rho)

JeV=1.6*1e-19;
me=9.11*1e-31;
c=3*1e8;
Erest=me*c^2/JeV;
tau=E./Erest;
% [beta,gamma,v]=vrel(E);
[beta,gamma]=betacalc(E);
v=beta.*c;

Na=6.02*1e23;
re=2.818*1e-15;
Ftau=1-beta.^2+(tau.^2/8-(2.*tau+1)*log(2))./(tau+1).^2;
delta=0;

const=2*pi*re^2*Erest*Na*Z./(beta.^2*Molmass);
Sp=rho*const.*(log(tau.^2.*(tau+2)/(2*(I/Erest)^2))+Ftau-delta);
Sp=Sp./1e9;