%% testbench for phonon scattering parameters
clear;
clc;
close all;

eps0=1.5^2;
epsinf=1;
E_optphonon=0.1;

E=0.2:0.01:10;

[imfp_T,imfp_creation,imfp_annihilation,theta]=ephscatt(E,eps0,epsinf,E_optphonon);

figure;
plot(E,imfp_T,'b','linewidth',3.0);

hold on;
plot(E,imfp_creation,'k','linewidth',3.0);

plot(E,imfp_annihilation,'r','linewidth',3.0);

set(gca,'XScale','linear','YScale','log');

%% testbench for phonon scattering parameters: angles
clear;
clc;
close all;

eps0=1.5^2;
epsinf=1;
hbarw=0.1;

E=1;
Eprime=E-hbarw;

x=hbarw/E;

ntrials=1000;
for i = 1:ntrials
    B=(E+Eprime+2*sqrt(E*Eprime))/(E+Eprime-2*sqrt(E*Eprime));
    rng('shuffle');
    randval=rand;
    cosval(i)=(2-x)/(2*sqrt(1-x)).*(1-B^randval)+B^randval;
end

theta=acos(cosval);

figure;
plot(theta.*180/pi);

figure;
hist(theta.*180/pi,100);

figure;
plot(cosval>1);
hold on;
plot(cosval<1,'k');
% theta=[0:0.1:90].*pi/180;
% 
% W=2*sqrt(E*Eprime)/log((E+Eprime+2*sqrt(E*Eprime))/(E+Eprime-2*sqrt(E*Eprime))).*sin(theta)./(E+Eprime-2*sqrt(E*Eprime).*cos(theta));
% 
% figure;
% plot(theta.*180/pi,W);

