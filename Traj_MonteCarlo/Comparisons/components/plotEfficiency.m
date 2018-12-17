function handleArray = plotEfficiency(resultObject,configuration)
    FIGURE_TEMP = 3300;
    %% Initialization
    incEngy         =   configuration.incEnergy;
    runsToCompare   =   configuration.runsToCompare;
    nTrials         =   configuration.nTrials;    
    nRuns = length(resultObject);   
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
    plot(energy,nAcids/2,'-x');
    title('Number of acids per electron');
    xlabel('Incident energy (eV)')
    
    figure(8001) 
    hold on
    plot(energy,nAcids./energy/2,'-x');    
    title('Number of acids per eV per electron');
    xlabel('Incident energy (eV)')
end