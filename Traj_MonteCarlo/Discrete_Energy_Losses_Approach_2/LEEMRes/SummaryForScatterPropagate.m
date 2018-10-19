%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is an anlysis script that summrizes the output of
%%% the scattering -> propagation framework
%%%                                                             --jonathan 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
set(0,'DefaultFigureWindowStyle','docked')
%% Loading the files (Suchit): 
counter = 1;
for i = 1:1000
    load(strcat('Center_thetaSuchit\Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
    %close(1);
    if size(acid_xyz,1)>=1
        spread_80(i)= sqrt(sum(acid_xyz(:,1).^2 +acid_xyz(:,2).^2)./size(acid_xyz,1));
        acid_xyz_80(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
counter = 1;
for i = 1:1000
    load(strcat('Center_thetaSuchit\Ein=30.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
    %close(1)
    if size(acid_xyz,1) >=1
        spread_30(i)= sqrt(sum(acid_xyz(:,1).^2 +acid_xyz(:,2).^2)./size(acid_xyz,1));
        acid_xyz_30(counter:counter + size(acid_xyz,1)-1,:)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end

%% Loading the files (jonathan): 
counter = 1;
for i = 1:1000
    load(strcat('Center_thetaJonathan\Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
    %close(1);
    if size(acid_xyz,1)>=1
        acid_xyz_80_J(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end


counter = 1;
for i = 1:1000
    load(strcat('Center_thetaJonathan\Ein=30.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
    %close(1)
    if size(acid_xyz,1) >=1
        acid_xyz_30_J(counter:counter + size(acid_xyz,1)-1,:)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end

%% Loading files produced after fixing the angular randomization in the differential X-sect
counter = 1;
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\LowEnergyRand\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_2(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
counter = 1;
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\LowEnergyRandReverted\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_3(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
counter = 1;
acid_xyz_80_J_4=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\LowEnergyRand_ArrayDecla_Reverted\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_4(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_2 = std(acid_xyz_80_J_2);
sig_80_J_3 = std(acid_xyz_80_J_3);
sig_80_J_4 = std(acid_xyz_80_J_4);
mu_80_J_2 = mean(acid_xyz_80_J_2,1);
mu_80_J_3 = mean(acid_xyz_80_J_3,1);
mu_80_J_4 = mean(acid_xyz_80_J_4,1);

counter = 1;
acid_xyz_80_J_5=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\LowERand_ArD_Reverted_Restarted\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_5(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_5 = std(acid_xyz_80_J_5);
mu_80_J_5 = mean(acid_xyz_80_J_5,1);

counter = 1;
acid_xyz_80_J_6=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\LowERand_ArD_Reverted_Restarted_2\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_6(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_6 = std(acid_xyz_80_J_6);
mu_80_J_6 = mean(acid_xyz_80_J_6,1);

%% Calculating the spread (the 2-D gaussian width)
com_80 = mean(acid_xyz_80,1);
sig_xy_80 = sqrt((sum((acid_xyz_80(:,1)-com_80(1)).^2)...
    +sum((acid_xyz_80(:,2)-com_80(2)).^2))/size(acid_xyz_80,1))
sig_xyz_80 = sqrt((sum((acid_xyz_80(:,1)-com_80(1)).^2)...
    +sum((acid_xyz_80(:,2)-com_80(2)).^2)...
    +sum((acid_xyz_80(:,3)-com_80(3)).^2))/size(acid_xyz_80,1))
com_30 = mean(acid_xyz_30,1);
sig_xy_30 = sqrt((sum((acid_xyz_30(:,1)-com_30(1)).^2)...
    +sum((acid_xyz_30(:,2)-com_30(2)).^2))/size(acid_xyz_30,1))
sig_xyz_30 = sqrt((sum((acid_xyz_30(:,1)-com_30(1)).^2)...
    +sum((acid_xyz_30(:,2)-com_30(2)).^2)...
    +sum((acid_xyz_30(:,3)-com_30(3)).^2))/size(acid_xyz_30,1))
%% Obtaining standard deviation in all directions
% Independent sigma for reference
sig_80 = std(acid_xyz_80)
sig_30 = std(acid_xyz_30)

sig_80_J = std(acid_xyz_80_J)
sig_30_J = std(acid_xyz_30_J)

%% Some graphs

r_80    =   sqrt((acid_xyz_80(:,1)-com_80(1)).^2+...
    (acid_xyz_80(:,2)-com_80(2)).^2+...
    (acid_xyz_80(:,3)-com_80(2)).^2);
r_30    =   sqrt((acid_xyz_30(:,1)-com_30(1)).^2+...
    (acid_xyz_30(:,2)-com_30(2)).^2+...
    (acid_xyz_30(:,3)-com_30(2)).^2);


% Histograms
edges = 0:0.5:10;
figure(401)
hold off
histogram(r_80,25,'Normalization','probability','BinEdges',edges)
hold on
histogram(r_30,20,'Normalization','probability','BinEdges',edges)
legend('80-eV','30-eV');
%% Histogram for the _J series
