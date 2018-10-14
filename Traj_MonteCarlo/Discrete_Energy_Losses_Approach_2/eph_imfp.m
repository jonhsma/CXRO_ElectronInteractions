function imfp_eph=eph_imfp(parms)
% imfp_T,imfp_creation,imfp_annihilation,theta % outputs from eph_imfp

% eps0=1.5^2;
% epsinf=1;
% E_optphonon=0.1; % eV;

E=parms.E;
eps0=parms.eps0;
epsinf=parms.epsinf;
E_optphonon=parms.hbarw;

imfp_eph=ephscatt(E,eps0,epsinf,E_optphonon);
imfp_eph=(imfp_eph).*1e9; % nm;
