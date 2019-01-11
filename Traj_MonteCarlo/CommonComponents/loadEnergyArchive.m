%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function loads the database and stiches the chunks together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stichedEnergyArchive = loadEnergyArchive(basepath,filename)
    fileCounter = 1;
    stichedEnergyArchive = [];
    while(exist(strcat(basepath,filename,'_',num2str(fileCounter),'.mat'),'file'))
        currChunk = load(strcat(basepath,filename,'_',num2str(fileCounter)));
        stichedEnergyArchive =  [stichedEnergyArchive currChunk.archiveChunk];
        fileCounter = fileCounter +1;
    end
end