%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is a quick diagnostic script after each run of TrajSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pull up a certain workspace from the data folder
clear;
nameOfRun = 'AcidEventPosLogging_0';
load(strcat('..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\',...
        nameOfRun,'\',...
        'WorkSpace.mat'))
%% Pull up the acid distribution
counter = 1;
counter_f = 1;
acid_xyz_accul=[];
for i = 1:1000
    try
        load(strcat('..\..\..\..\',...
            'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\',...
            nameOfRun,...
            '\Ein=80.00_Dose=32.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
        if size(acid_xyz,1)>=1
            acid_xyz_accul(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
            counter = counter+size(acid_xyz,1);
        end
        if size(acid_fine_xyz,1)>=1
            acid_fine_xyz_accul(counter:counter + size(acid_fine_xyz,1)-1,1:3)=acid_fine_xyz(:,:);
            counter_f = counter_f+size(acid_fine_xyz,1);
        end
        if counter_f~=counter
            fprintf('Length mismatch between pixelated acid positions and the continuous one')
        end        
    catch exception
        fprintf('Problem reading the %d-th file\n', i);
        break;
    end
end

%% Acid Statistics
sig_acid        = std(acid_xyz_accul);
mu_acid         = mean(acid_xyz_accul,1);
mu_acid_fine    = mean(acid_fine_xyz_accul,1);

%%% Sptial distribution of acids
figure(5001);
hold off
histogram(acid_xyz_accul(:,1)-mu_acid(1),21,'BinEdges',-5.25:0.5:5.25);
hold on
histogram(acid_xyz_accul(:,2)-mu_acid(2),21,'BinEdges',-5.25:0.5:5.25);
histogram(acid_xyz_accul(:,3)-mu_acid(3),21,'BinEdges',-5.25:0.5:5.25);
title('Distribution of acid positions');
legend('x','y','z')

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


figure(5004);
hold off
histogram(r3,'BinEdges',-0.125:0.25:5.125);
hold on
histogram(r3_fine,'BinEdges',-0.0625:0.125:5.0625);
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
figure(4004);hist(theta_global,binsX);
title('Theta distribution for all events in scattering frame');
%% Electron Statistics
fprintf('Standard deviations and means of data\n');
sampleSize = floor(size(xyz_electron_global,1)/3)*3;
disp(std(xyz_electron_global));
disp(mean(xyz_electron_global));
fprintf('Standard deviations of subsets \n');
disp(std(xyz_electron_global(1:sampleSize/3,:)))
disp(std(xyz_electron_global(sampleSize/3+1:sampleSize/3*2,:)))
disp(std(xyz_electron_global(sampleSize/3*2+1:sampleSize,:)))

%%