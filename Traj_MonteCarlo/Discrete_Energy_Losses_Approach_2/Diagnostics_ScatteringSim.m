%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is a quick diagnostic script after each run of TrajSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pull up a certain workspace from the data folder
clear;
nameOfRun = 'NoCoarseGrain_2_OffCntr';
load(strcat('..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\',...
        nameOfRun,'\',...
        'WorkSpace.mat'));
res     =   0.25;
doseStr =   '0.00';
incEngy =   '80.00';
trialN =   1000;
set(0,'DefaultFigureWindowStyle','docked')
%% Pull up the acid distribution
counter = 1;
counter_f = 1;
acid_xyz_accul=[];
acid_fine_xyz_accul=[];
for i = 1:trialN
    try
        load(strcat('..\..\..\..\',...
            'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\',...
            nameOfRun,...
            '\Ein=',incEngy,'_Dose=',doseStr,...
            'epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
        if size(acid_xyz,1)>=1
            acid_xyz_accul(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
            counter = counter+size(acid_xyz,1);
        end
        if exist('acid_fine_xyz','var')==1 && size(acid_fine_xyz,1)>=1
            acid_fine_xyz_accul(counter_f:counter_f + size(acid_fine_xyz,1)-1,1:3)=acid_fine_xyz(:,:);
            counter_f = counter_f+size(acid_fine_xyz,1);
        end
        if counter_f~=counter && exist('acid_fine_xyz','var')==1
            fprintf('Length mismatch between pixelated acid positions and the continuous one\n')
        end        
    catch exception
        disp(exception);
        fprintf('Problem reading the %d-th file\n', i);
        break;
    end
end
%% Acid Statistics
sig_acid        = std(acid_xyz_accul);
mu_acid         = mean(acid_xyz_accul,1);
mu_acid_fine    = mean(acid_fine_xyz_accul,1);

sampleSize = floor(size(acid_xyz_accul,1)/3)*3;
fprintf('Standard Deviations and means of acid positions\n');
disp(sig_acid)
disp(mu_acid)
fprintf('Standard deviations of subsets \n');
disp(std(acid_xyz_accul(1:sampleSize/3,:)))
disp(std(acid_xyz_accul(sampleSize/3+1:sampleSize/3*2,:)))
disp(std(acid_xyz_accul(sampleSize/3*2+1:sampleSize,:)))


%%% Sptial distribution of acids
binEdges = -5-res/2:res:5+res/2;

figure(5001);
subplot(3,1,1)
hold off
histogram(acid_xyz_accul(:,1)-mu_acid(1),'BinEdges',binEdges);
title('Distribution of acid positions-pixelated coorinates');
legend('x')
subplot(3,1,2)
histogram(acid_xyz_accul(:,2)-mu_acid(2),'BinEdges',binEdges);
legend('y')
subplot(3,1,3)
histogram(acid_xyz_accul(:,3)-mu_acid(3),'BinEdges',binEdges);
legend('z')

%%% Angle distribution of acids
r3      = acid_xyz_accul(:,:)-mu_acid(:)';
r3_fine = acid_fine_xyz_accul(:,:)-mu_acid_fine(:)';


rabs        =   sqrt(r3(:,1).^2+r3(:,2).^2+r3(:,3).^2);
rabs_fine   =   sqrt(r3_fine(:,1).^2+r3_fine(:,2).^2+r3_fine(:,3).^2);


acid_direction = r3./rabs;
acid_fine_direction = r3_fine./rabs_fine;

figure(5002)
hold off
scatter3(acid_direction(:,1),...
    acid_direction(:,2),...
    acid_direction(:,3),'.')
title('Distribution of direction of acids-pixelated coordinates');
pbaspect([1 1 1]);

figure(5003)
hold off
scatter3(acid_fine_direction(:,1),...
    acid_fine_direction(:,2),...
    acid_fine_direction(:,3),'.')
title('Distribution of direction of acids-machine precision');
pbaspect([1 1 1]);

radialBins = 0:res/2:5;

figure(5004);
hold off
histogram(rabs,'BinEdges',radialBins);
hold on
histogram(rabs_fine,'BinEdges',radialBins);
title('Radial Distribution of acids');
legend('pixelated acid position','machine precision event position')
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
figure(4004);hist(theta_global,-pi:pi/8:pi);
title('Theta distribution for all events in scattering frame');
%% Electron Statistics
fprintf('Standard deviations and means of electron event positions\n');
sampleSize = floor(size(xyz_electron_global,1)/3)*3;
disp(std(xyz_electron_global));
disp(mean(xyz_electron_global));
fprintf('Standard deviations of subsets \n');
disp(std(xyz_electron_global(1:sampleSize/3,:)))
disp(std(xyz_electron_global(sampleSize/3+1:sampleSize/3*2,:)))
disp(std(xyz_electron_global(sampleSize/3*2+1:sampleSize,:)))

mu_e = mean(xyz_electron_global);

figure(4005);
hold off
histogram(xyz_electron_global(:,1)-mu_e(1),21,'BinEdges',-5.25:0.5:5.25);
hold on
histogram(xyz_electron_global(:,2)-mu_e(2),21,'BinEdges',-5.25:0.5:5.25);
histogram(xyz_electron_global(:,3)-mu_e(3),21,'BinEdges',-5.25:0.5:5.25);
title('Distribution of electron scattering event positions');
legend('x','y','z')
%% Event Activation Analysis
diff3   = acid_fine_xyz_accul-acid_xyz_accul;
diffabs     =   sqrt(diff3(:,1).^2+diff3(:,2).^2+diff3(:,3).^2);
figure(3000); 
plot(rabs_fine,rabs,'.')
title('Event and activation radial positions')
xlabel('Event radial position')
ylabel('Activated acid voxel center radial position')
axis([0,5,0,5])

figure(3001)
plot(rabs_fine,diffabs,'.')
title('Distance of activation')
xlabel('Event radial position')
ylabel('Distance between acid voxel center and activation event')
axis([0,3,0,2.5])