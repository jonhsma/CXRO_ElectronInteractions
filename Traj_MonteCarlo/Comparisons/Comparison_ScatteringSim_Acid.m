%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script compares multiple runs with the scatter sim
%%% Acid stats and only acid stats are used because loading an entire
%%% workspace carries the risk of overwriting existing variables
%%% --> On the todo list maybe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resultObject=Comparison_ScatteringSim_Acid(varargin)
    addpath('components\')
    %% Select the runs to compare
    switch nargin
        case 0
            %%% Defaults in case the code is run without arguments
            runsToCompare =...
                {'NoCoarseGrain_0_Calib';
                'Ref_20181103'%'RealScatt_80_2_1nmGrid'};
                };
            nRuns = size(runsToCompare,1);
            incEngy =...
                {'80.00';'80.00'};
            doseStr =...
                {'0.00';'0.00'};
            nTrials = {1000,1000};
        case 1
            configuration   =   varargin{1};
            runsToCompare   =   configuration.runsToCompare;
            incEngy         =   configuration.incEnergy;
            doseStr         =   configuration.doseStr;
            nTrials         =   configuration.nTrials;
            nRuns = size(runsToCompare,1);
    end
   
    %% Pull up the acid distribution
    resultObject = LoadArchives(runsToCompare,nRuns,incEngy,nTrials,doseStr);
    %% Display Settings
    %%% Setting the resolution of the histograms
    %%% Run the defaults first
    res = 0.25;    
    cuOption    =   'cdf';
    distOption  =   'count';
    legendNotes =   runsToCompare;
    rangeLimit  =   5;
    if exist('configuration','var')
        if isfield(configuration,'res')
            res         =   configuration.res;
        end
        if isfield(configuration,'cuOption')
            cuOption    =   configuration.cuOption;
        end
        if isfield(configuration,'distOption')
            distOption  =   configuration.distOption;
        end
        if isfield(configuration,'rangeLimit')
            rangeLimit  =   configuration.rangeLimit;
        end
    end
    
    for ii = 1:nRuns
        target             = strrep(runsToCompare{ii},'_','\_');
        legendNotes{ii} = strcat(target,'; N_{trials}=',num2str(nTrials{ii}));
    end
    %% Acid Statistics
    plotAcidDistCurves(resultObject,configuration)
    %{
    %% Acid linear Histograms
    binEdges = -rangeLimit-res/2:res:rangeLimit+res/2;
    figure(5010);

    for iDim = 1:3
        subplot(3,1,iDim)
        hold off
        for ii=1:nRuns
            this = resultObject{ii};
            histogram(this.acid_xyz_accul(:,iDim),...
                'BinEdges',binEdges,'Normalization',distOption);
            hold on
        end
        if iDim ==1
            title('Distribution of acid positions-pixelated coorinates');
            legend(legendNotes)
        end
    end
    xlabel('position (nm)')
    %% Acid radial histograms
    binEdges = 0:res:2*rangeLimit;
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
            'BinEdges',binEdges,'Normalization',cuOption);
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
            'BinEdges',binEdges,'Normalization',cuOption);
        hold on
    end
    title('Cumulatice Radial Distribution of Activation Events')
    xlabel('Distance from Origin (nm)')
    legend(legendNotes,'Location','southeast')
    %}
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
    axis([0,rangeLimit,0,rangeLimit])

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
    axis([0,rangeLimit,0,rangeLimit/4])
    %% Tabulation
    tabulateSummary(resultObject);
end