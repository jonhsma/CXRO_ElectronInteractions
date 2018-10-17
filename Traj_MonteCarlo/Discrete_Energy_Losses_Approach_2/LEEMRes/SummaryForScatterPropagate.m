%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is an anlysis script that summrizes the output of
%%% the scattering -> propagation framework
%%%                                                             --jonathan 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading the files
clear;
for i = 1:5
    load(strcat('Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
    spread(i)= sqrt(sum(acid_xyz(:,1).^2 +acid_xyz(:,2).^2)./size(acid_xyz,1));
end
