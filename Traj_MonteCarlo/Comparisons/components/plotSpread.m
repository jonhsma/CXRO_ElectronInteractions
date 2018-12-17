function handleArray = plotSpread(resultObject,configuration)
    FIGURE_TEMP = 3300;
    %% Initialization
    incEngy         =   configuration.incEnergy; 
    nRuns = length(resultObject);   
    %% Getting the numbers
    energy = zeros([1 nRuns]);
    rRMS = energy;
    
    for ii = 1:nRuns
       energy(ii)   =   str2double(incEngy{ii});
       rRMS(ii)   =   sqrt(mean(sum(resultObject{ii}.acid_xyz_accul.^2,2)));
       rEXP(ii)   =   mean(sqrt(sum(resultObject{ii}.acid_xyz_accul.^2,2)));
    end    
    %% Plots    
    figure(9000) 
    hold on
    plot(energy,rRMS);
    title('<r^2>^{1/2}');
    xlabel('Incident energy (eV)')
    figure(9001) 
    hold on
    plot(energy,rEXP);
    title('<r>');
    xlabel('Incident energy (eV)')
end