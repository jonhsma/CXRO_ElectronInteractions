%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script test if modularization codes are working as expected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MFP generation
addpath('..\');
optdata_path='..\..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';
scattdata.optical=load(['..\' optdata_path 'Sp_Fuji_IMFP_Inelastic_Components_Ef=15.5eV_Elossmin=0.001eV_Erange=[16,200]_DDCSData.mat']);
pathname='..\..\Traj_MonteCarlo\Discrete_Energy_Losses_Approach_2\DDCSData\';
%DCS data
filename='DDCSdata_Fuji_Ef=15.5_Elossmin=0.001eV_Erange=[16,200].mat';
scattdata.optical.inel_dcsdata  = load(['..\' pathname filename]);
scattdata.E_inel_thr            = min(scattdata.optical.E);

e_test = 20:0.25:100;
counter = 0;
MFP_original    = zeros([1 length(e_test)]);
MFP_modular      = MFP_original; 

event.lowEthr   =   20;
event.lowEimfp  =   3.67;

for e_i = e_test
        counter = counter+1;
    
        controlparm.onlyimfp=1;
        if isfield(scattdata.optical,'inel_dcsdata')
            Elossrand_opt=genrandEloss_OptData_JHM(scattdata.optical,e_i,scattdata.optical.inel_dcsdata,controlparm);
        else
            Elossrand_opt=genrandEloss_OptData_JHM(scattdata.optical,e_i,controlparm);
        end
        if e_i>event.lowEthr
            controlparm.onlyimfp=0;
            imfp_opt=Elossrand_opt.imfp; % imfp in nm; the above lines only calculate imfp due to the controlparms.onlyimfp line abvove.
            if isnan(imfp_opt)
                imfp_opt=Inf;
            end
        else
            imfp_opt=event.lowEimfp;
        end
    MFP_original(counter)    =   imfp_opt; 
   
    MFP_modular(counter) = genMFP_OptData(event,scattdata.optical,e_i);
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