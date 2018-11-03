%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script compares multiple runs with the scatter sim
%%% Acid stats and only acid stats are used because loading an entire
%%% workspace carries the risk of overwriting existing variables
%%% --> On the todo list maybe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select the runs to compare
runsToCompare =...
    {'NoCoarseGrain_0_Calib'};...
    %'RealScatt_80_1_1nmGrid''RealScatt_80_2_1nmGrid'};
    %};
nRuns = size(runsToCompare,1);
incEngy =...
    {'80.00'};%;'80.00'};
doseStr =...
    {'0.00'};%;'2.00'};
resultObject = cell([1 nRuns]);


%% Display Settings
%%% Setting the resolution of the histograms
res = 0.25;
legendNotes = runsToCompare;
for ii = 1:nRuns
    target             = strrep(runsToCompare{ii},'_','\_');
    legendNotes{ii} = target;
end
cumOption   =   'cdf';
distOption  =   'count';
%% Pull up the acid distribution
for runCnt = 1:size(runsToCompare,1)
    counter = 1;
    counter_f = 1;
    resultObject{runCnt}.acid_xyz_accul=[];
    resultObject{runCnt}.acid_fine_xyz_accul=[];
    
    fprintf('Currently Processing')
    disp(runsToCompare(runCnt,:));
    
    for i = 1:1000
        try
            load(strcat('..\..\..\..\',...
                'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\',...
                runsToCompare{runCnt},...
                '\Ein=',incEngy{runCnt},...
                '_Dose=',doseStr{runCnt},...
                'epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
            if size(acid_xyz,1)>=1
                resultObject{runCnt}.acid_xyz_accul(counter:counter + size(acid_xyz,1)-1,1:3)...
                    =acid_xyz(:,:);
                counter = counter+size(acid_xyz,1);
            end
            if size(acid_fine_xyz,1)>=1
                resultObject{runCnt}.acid_fine_xyz_accul(counter_f:counter_f + size(acid_fine_xyz,1)-1,1:3)...
                    =acid_fine_xyz(:,:);
                counter_f = counter_f+size(acid_fine_xyz,1);
            end
            if counter_f~=counter
                fprintf('Length mismatch between pixelated acid positions and the continuous one')
            end        
        catch exception
            disp(exception);
            fprintf('Problem reading the %d-th file\n', i);
            break;
        end
    end
    
    %%% The means
    resultObject{runCnt}.mu_acid =...
        mean(resultObject{runCnt}.acid_xyz_accul);
    resultObject{runCnt}.mu_acid_fine =...
        mean(resultObject{runCnt}.acid_fine_xyz_accul);
    fprintf('The mean of position of acids\n');
    disp(resultObject{runCnt}.mu_acid);
    
    %%% The spreads
    resultObject{runCnt}.sig_acid =...
        std(resultObject{runCnt}.acid_xyz_accul);
    resultObject{runCnt}.sig_acid_fine =...
        std(resultObject{runCnt}.acid_fine_xyz_accul);
    fprintf('The standard deviation in position of acids\n');
    disp(resultObject{runCnt}.sig_acid);
    fprintf('The standard deviation in position of activation events\n');
    disp(resultObject{runCnt}.sig_acid_fine);    
    
    %%% Corrected positions
    resultObject{runCnt}.r3=...
        resultObject{runCnt}.acid_xyz_accul...
        -resultObject{runCnt}.mu_acid;
    resultObject{runCnt}.r3_fine=...
        resultObject{runCnt}.acid_fine_xyz_accul...
        -resultObject{runCnt}.mu_acid_fine;
    
    %%% The radii
    resultObject{runCnt}.rabs =...
        sqrt(sum(resultObject{runCnt}.r3.^2,2));
    resultObject{runCnt}.rabs_fine =...
        sqrt(sum(resultObject{runCnt}.r3_fine.^2,2));
    
    %%% The event-activation distances
    resultObject{runCnt}.acidEventdist = ...
        sqrt(sum((resultObject{runCnt}.r3_fine-resultObject{runCnt}.r3).^2,2));
end
%% Acid linear Histograms
binEdges = -5-res/2:res:5+res/2;
figure(5010);

for iDim = 1:3
    subplot(3,1,iDim)
    hold off
    for ii=1:nRuns
        this = resultObject{ii};
        histogram(this.acid_xyz_accul(:,iDim),...
            'BinEdges',binEdges,'Normalization','probability');
        hold on
    end
    if iDim ==1
        title('Distribution of acid positions-pixelated coorinates');
        legend(legendNotes)
    end
end
xlabel('position (nm)')
%% Acid radial histograms
binEdges = 0:res:10;
figure(5011);
hold off
for ii=1:nRuns
    this = resultObject{1,ii};
    histogram(this.rabs(:,1),...
        'BinEdges',binEdges,'Normalization',distOption);
    hold on
end
title('RadialDistribution of acids')
xlabel('Radial position (nm)')
legend(legendNotes,'Location','northeast')

figure(5012);
hold off
for ii=1:nRuns
    this = resultObject{1,ii};
    histogram(this.rabs(:,1),...
        'BinEdges',binEdges,'Normalization',cumOption);
    hold on
end
title('Culmulative Radial Counts of Acids')
xlabel('Radial position (nm)')
legend(legendNotes,'Location','northwest')

figure(5013);
hold off
for ii=1:nRuns
    this = resultObject{1,ii};
    histogram(this.rabs_fine(:,1),...
        'BinEdges',binEdges,'Normalization',distOption);
    hold on
end
title('RadialDistribution of Activation Events')
xlabel('Distance from Origin (nm)')
legend(legendNotes,'Location','northeast')



figure(5014);
hold off
for ii=1:nRuns
    this = resultObject{1,ii};
    histogram(this.rabs_fine(:,1),...
        'BinEdges',binEdges,'Normalization',cumOption);
    hold on
end
title('Cumulatice Radial Distribution of Activation Events')
xlabel('Distance from Origin (nm)')
legend(legendNotes,'Location','southeast')
%% Activation Statistics
figure(3010);
hold off
for ii = 1:nRuns
    this = resultObject{ii};
    plot(this.rabs_fine,this.rabs,'.')
    hold on
end
title('Event and activation radial positions')
xlabel('Event distance from origin')
ylabel('Distance from origin - activated acid voxel center')
legend(legendNotes,'Location','northwest')
axis([0,5,0,5])

figure(3011)
hold off
for ii = 1:nRuns
    this = resultObject{ii};
    plot(this.rabs_fine,this.acidEventdist,'.')
    hold on;
end
title('Distance of activation')
xlabel('Event distance from orgin')
ylabel('Distance between acid voxel center and activation event')
axis([0,3,0,2.5])