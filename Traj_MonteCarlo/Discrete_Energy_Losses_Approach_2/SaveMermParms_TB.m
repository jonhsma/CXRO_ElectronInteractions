%% testbench for saving the Mermin parameters
clear;
clc;
close all;

help SaveMermParms

dest_dir='MerminParameters\';
fname=[dest_dir 'MermParms_Fuji1201E_EELS.mat'];

physparms.rho=1.19;
physparms.Mw=120;
physparms.Ef=15.5;
physparms.q=1e9;
physparms.hbar=6.626*1e-34/(2*pi);

mermparms_lsqfit=[];
mermparms_lsqfit=[mermparms_lsqfit;[0.051 6.22 3.14]];
mermparms_lsqfit=[mermparms_lsqfit;[0.295 23.22 10.69]];
mermparms_lsqfit=[mermparms_lsqfit;[0.224 32.58 26.49]];
mermparms_lsqfit=[mermparms_lsqfit;[0.001 295.7 24]];

mermparms_legend='[Ai Ei gamma_i]';

save(fname,'physparms','mermparms_lsqfit','mermparms_legend');