%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the distributions of a result object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handleArray = plotAcidDistCurves(resultObject,varargin)
    FIGURE_TEMP = 3300;
    %% Initialization
    if nargin == 2
        configuration = varargin{1};
        incEngy         =   configuration.incEnergy;
        runsToCompare   =   configuration.runsToCompare;
        nTrials         =   configuration.nTrials;
    end
    
    nRuns = length(resultObject);
    %% Display Settings
    %%% Setting the resolution of the histograms
    %%% Run the defaults first
    res = 0.25;    
    cuOption    =   'cdf';
    distOption  =   'count';
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
    
    if exist('runsToCompare','var')
        for ii = 1:nRuns
            target             = strrep(runsToCompare{ii},'_','\_');
            legendNotes{ii} = strcat(target(1:8),...
                ';',incEngy{ii},...
                '; N_{trials}=',num2str(nTrials{ii}));
        end
    else
           for ii = 1:nRuns           
            legendNotes{ii} = strcat('Item-',num2str(ii));
        end
    end
    %% Acid linear Histograms
    binEdges = -rangeLimit-res/2:res:rangeLimit+res/2;
    binCenter= -rangeLimit:res:rangeLimit;
    figure(5010);

    for iDim = 1:3
        subplot(3,1,iDim)
        hold off
        for ii=1:nRuns
            this = resultObject{ii};
            figure(FIGURE_TEMP)
            h = histogram(this.acid_xyz_accul(:,iDim),...
                'BinEdges',binEdges,'Normalization',distOption);
            figure(5010)
            plot(binCenter,h.Values);
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
    binCenter = res/2:res:2*rangeLimit-res/2;
    figure(5011);
    hold off
    for ii=1:nRuns
        this = resultObject{1,ii};
        figure(FIGURE_TEMP);
        h = histogram(this.rabs(:,1),...
            'BinEdges',binEdges,'Normalization',distOption);
        figure(5011)
        plot(binCenter,h.Values);
        hold on
    end
    title('RadialDistribution of acids')
    xlabel('Radial position (nm)')
    legend(legendNotes,'Location','northeast')

    figure(5012);
    hold off
    for ii=1:nRuns
        this = resultObject{1,ii};
        figure(FIGURE_TEMP)
        h = histogram(this.rabs(:,1),...
            'BinEdges',binEdges,'Normalization',cuOption);
        figure(5012)
        plot(binCenter,h.Values);
        hold on
    end
    title('Culmulative Radial Counts of Acids')
    xlabel('Radial position (nm)')
    legend(legendNotes,'Location','northwest')

    figure(5013);
    hold off
    for ii=1:nRuns
        this = resultObject{1,ii};
        figure(FIGURE_TEMP)
        h = histogram(this.rabs_fine(:,1),...
            'BinEdges',binEdges,'Normalization',distOption);
        figure(5013)
        plot(binCenter,h.Values);
        hold on
    end
    title('RadialDistribution of Activation Events')
    xlabel('Distance from Origin (nm)')
    legend(legendNotes,'Location','northeast')

    figure(5014);
    hold off
    for ii=1:nRuns
        this = resultObject{1,ii};
        figure(FIGURE_TEMP);
        h = histogram(this.rabs_fine(:,1),...
            'BinEdges',binEdges,'Normalization',cuOption);
        figure(5014)
        plot(binCenter,h.Values);
        hold on
    end
    title('Cumulatice Radial Distribution of Activation Events')
    xlabel('Distance from Origin (nm)')
    legend(legendNotes,'Location','southeast')
end