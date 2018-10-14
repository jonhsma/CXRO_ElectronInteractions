function [beta,gamma]=betacalc(ke)
%%% [beta, gamma] = betaclac(ke);
%%% ke: kinetic energy in eV

JeV=1.6*1e-19;
mo=9.11*1e-31;
c=3*1e8;
relmass=(ke.*JeV+mo*c^2)./c^2;

gamma=relmass./mo;
beta=sqrt(1-1./gamma.^2);