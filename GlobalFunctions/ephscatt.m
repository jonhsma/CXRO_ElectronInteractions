function [imfp_T,imfp_creation,imfp_annihilation,theta]=ephscatt(E,eps0,epsinf,E_optphonon)

hbarw=E_optphonon; % eV, used by Dapor 2012
kT=0.025; % eV, assuming room temp.
ao=53*1e-12; % SI units

nT=1./(exp(hbarw/kT)-1);

const=(nT+1)*(1/epsinf-1/eps0)/(2*ao);
num=1+sqrt(1-hbarw./E);
den=1-sqrt(1-hbarw./E);
laminv=const*hbarw./E.*log(num./den);

const2=(nT)*(1/epsinf-1/eps0)/(2*ao);
num2=1+sqrt(1+hbarw./E);
den2=-1+sqrt(1+hbarw./E);
laminv2=const2*hbarw./E.*log(num2./den2);

imfp_T=1e9./(laminv+laminv2); % convert to nm with 1e9
imfp_creation=1e9./laminv; % convert to nm with 1e9
imfp_annihilation=1e9./laminv2; % convert to nm with 1e9

%%% angle analytical from Villarubia's result

%%%% THESE 4 LINES worked correctly [if you need it in the future]
% B=(E+Eprime+2*sqrt(E*Eprime))/(E+Eprime-2*sqrt(E*Eprime));
% rng('shuffle');
% randval=rand;
% cosval(i)=(2-x)/(2*sqrt(1-x)).*(1-B^randval)+B^randval;
%%%% End of THE 4 LINES worked correctly [if you need it in the future]

theta=[];
% for i = 1:length(E)
%     x=hbarw./E(i);
    x=hbarw./E;
    B=(2-x+2.*sqrt(1-x))./(2-x-2.*sqrt(1-x));
    rng('shuffle');
    randval=rand;
    cosval=(2-x)./(2.*sqrt(1-x)).*(1-B.^randval) + B.^randval;
% %     if cosval>1
% %         dbg=1;
% %     end
%     theta(i)=acos(cosval);
    theta=acos(cosval);
% end

% imfp=1./(laminv+laminv2);
% flag=1;