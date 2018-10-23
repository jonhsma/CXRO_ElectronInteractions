%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is a quick diagnostic script after each run of TrajSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pull up a certain workspace from the data folder
clear;
load(strcat('..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\Quarter_nm_pixel_0\',...
        'WorkSpace.mat'))
%% Distributions of angles
figure(4000); hist(theta_init,20)
title('Incident electrons theta distribution');
figure(4001);scatter3(scattVector(1,:),scattVector(2,:),scattVector(3,:),'.')
title('post-scattering electron direction distribution');
daspect([1 1 1])
figure(4002);hist(secSpawningTheta,20)
title('Secondary electron sapwning distribution');
figure(4003);hist(thetaLog,20)
title('theta distribution for all events');
nBin = 7;
binsX = round((pi/nBin/2:pi/nBin:pi-pi/nBin/2)*10)/10;
figure(4004);hist(theta_global,binsX);
title('theta distribution for all events in scattering frame');
%% Basic Statistics
fprintf('Standard deviations and means of data\n');
sampleSize = floor(size(xyz_electron_global,1)/3)*3;
disp(std(xyz_electron_global));
disp(mean(xyz_electron_global));
fprintf('Standard deviations of subsets \n');
disp(std(xyz_electron_global(1:sampleSize/3,:)))
disp(std(xyz_electron_global(sampleSize/3+1:sampleSize/3*2,:)))
disp(std(xyz_electron_global(sampleSize/3*2+1:sampleSize,:)))

%%