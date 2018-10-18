%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is an anlysis script that summrizes the output of
%%% the scattering -> propagation framework
%%%                                                             --jonathan 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading the files (Suchit): 
clear;
counter = 1;
set(0,'DefaultFigureWindowStyle','docked')

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
