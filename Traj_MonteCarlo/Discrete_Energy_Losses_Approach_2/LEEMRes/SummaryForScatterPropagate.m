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

%% Isotropic random scattering (in lab frame)

counter = 1;
acid_xyz_80_J_7=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\IsoScatt_NoAngleAddition_0\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_7(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end

sig_80_J_7 = std(acid_xyz_80_J_7);
mu_80_J_7 = mean(acid_xyz_80_J_7,1);

counter = 1;
acid_xyz_80_J_8=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\IsoScatt_NoAngleAddition_1\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_8(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_8 = std(acid_xyz_80_J_8);
mu_80_J_8 = mean(acid_xyz_80_J_8,1);

%% Coordinate Transform Debug
counter = 1;
acid_xyz_80_J_9=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\RealScatt_CoordTransform_0\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_9(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_9 = std(acid_xyz_80_J_9);
mu_80_J_9 = mean(acid_xyz_80_J_9,1);

counter = 1;
acid_xyz_80_J_10=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\RealScatt_CoordTransform_1\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_10(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_10 = std(acid_xyz_80_J_10);
mu_80_J_10 = mean(acid_xyz_80_J_10,1);

counter = 1;
acid_xyz_80_J_11=[];
for i = 1:300
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\RealScatt_CoordTransform_2\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_11(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_11 = std(acid_xyz_80_J_11);
mu_80_J_11 = mean(acid_xyz_80_J_11,1);

counter = 1;
acid_xyz_80_J_12=[];
for i = 1:2000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\RealScatt_CoordTransform_3\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_12(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_12 = std(acid_xyz_80_J_12);
mu_80_J_12 = mean(acid_xyz_80_J_12,1);

%% Using Forward scattering to further narrow down the problem
%%% Later on an artangent problem was identified
counter = 1;
acid_xyz_80_J_13=[];
for i = 1:2000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\ForwardScatt_CoordTransform_0\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_13(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_13 = std(acid_xyz_80_J_13);
mu_80_J_13 = mean(acid_xyz_80_J_13,1);

counter = 1;
acid_xyz_80_J_15=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\ForwardScatt_CoordTransform_d2\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_15(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_15 = std(acid_xyz_80_J_15);
mu_80_J_15 = mean(acid_xyz_80_J_15,1);

%% Forward Scattering with all known mathematical problem fixed
counter = 1;
acid_xyz_80_J_16=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\ForwardScatt_CoordTransform_1\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_16(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_16 = std(acid_xyz_80_J_16);
mu_80_J_16 = mean(acid_xyz_80_J_16,1);

%% Orthogonal Scattering for comparison with forward scattering
counter = 1;
acid_xyz_80_J_17=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\OrthogonalScatt_0\',...
        'Ein=80.00_Dose=2.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_17(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_17 = std(acid_xyz_80_J_17);
mu_80_J_17 = mean(acid_xyz_80_J_17,1);
%% Orthogonal Scattering for comparison with forward scattering
counter = 1;
acid_xyz_80_J_18=[];
for i = 1:1000
    load(strcat('..\..\..\..\..\',...
        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\Quarter_nm_pixel_0\',...
        'Ein=80.00_Dose=32.00epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
   if size(acid_xyz,1)>=1
        acid_xyz_80_J_18(counter:counter + size(acid_xyz,1)-1,1:3)=acid_xyz(:,:);
        counter = counter+size(acid_xyz,1);
    end
end
sig_80_J_18 = std(acid_xyz_80_J_18);
mu_80_J_18 = mean(acid_xyz_80_J_18,1);

figure(5001);
hold off
histogram(acid_xyz_80_J_18(:,1)-mu_80_J_18(1),21,'BinEdges',-5.25:0.5:5.25);
hold on
histogram(acid_xyz_80_J_18(:,2)-mu_80_J_18(2),21,'BinEdges',-5.25:0.5:5.25);
histogram(acid_xyz_80_J_18(:,3)-mu_80_J_18(3),21,'BinEdges',-5.25:0.5:5.25);
legend('x','y','z')

%%% Angle distribution of acids
r3_18 = acid_xyz_80_J_18(:,:)-mu_80_J_18(:)';
r_18  = sqrt(r3_18(:,1).^2+r3_18(:,2).^2+r3_18(:,3).^2);
acid_direction_18 = r3_18./r_18;
figure(5002)
hold off
scatter3(acid_direction_18(:,1),...
    acid_direction_18(:,2),...
    acid_direction_18(:,3),'.')

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
