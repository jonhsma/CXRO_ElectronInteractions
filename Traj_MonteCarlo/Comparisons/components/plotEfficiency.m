function handleArray = plotEfficiency(resultObject,configuration,varargin)
    FIGURE_TEMP = 3300;
    %% Initialization
    incEngy         =   configuration.incEnergy;
    runsToCompare   =   configuration.runsToCompare;
    nTrials         =   configuration.nTrials;    
    nRuns = length(resultObject);   
    if nargin >=3
        electronsPerTrial = varargin{1};
    else
        electronsPerTrial = 1;
    end
    %% Getting the numbers
    energy = zeros([1 nRuns]);
    nAcids = energy;
    
    for ii = 1:nRuns
       energy(ii)   =   str2double(incEngy{ii});
       nAcids(ii)   =   size(resultObject{ii}.acid_xyz_accul,1)/nTrials{ii};
    end    
    %% Plots    
    figure(8000) 
    hold on
    plot(energy,nAcids/electronsPerTrial,'-x');
    title('Number of acids per electron');
    xlabel('Incident energy (eV)')
    
    figure(8001) 
    hold on
    plot(energy,nAcids./energy/electronsPerTrial,'-x');    
    title('Number of acids per eV per electron');
    xlabel('Incident energy (eV)')
end