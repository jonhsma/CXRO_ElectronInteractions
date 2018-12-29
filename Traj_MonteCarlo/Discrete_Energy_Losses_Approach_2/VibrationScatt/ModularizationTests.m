%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script test if modularization codes are working as expected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MFP generation
addpath('..\');
scattdata.vibr=load(['..\' cross_sect_data_path 'VibrExcit_Data_Khakoo_2013.mat']);
moleculeDensity = 1.2/120*6.02*1e23*10;

e_test = 5:0.25:20;
counter = 0;
MFP_original    = zeros([1 length(e_test)]);
MFP_modular      = MFP_original;  
for e_i = e_test
    counter = counter + 1;
    result   =   genrandEloss_Vibr(scattdata.vibr,e_i);
    invImfp_invCM      =   result.ics*moleculeDensity;
    MFP_original(counter)    =   1/invImfp_invCM*1e7; % imfp in nm; 
   
    MFP_modular(counter) = genMFP_Vibr(scattdata.vibr,e_i,moleculeDensity);
end
figure(7210)
hold off
plot(MFP_original, MFP_modular,'Xk');
hold on    
plot(MFP_modular, MFP_modular,'-r');
title({'MFP (in nm) comparison between two routines.';' Sould be identical'});
xlabel('Modular generation')
ylabel('Original Approach')
legend('Data','1:1 reference');