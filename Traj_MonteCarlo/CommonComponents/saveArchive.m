%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function saves the energyScanArchive in chunks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = saveArchive(energyScanArchive,basepath,filename)
    nTrials = length(energyScanArchive);
    accuNevents = 0;
    startingIdx = 1;
    chunkIdx    =   0;
    
    for ii = 1 : nTrials
        nIncidences = length(energyScanArchive{ii}.incidences);
        for jj = 1:nIncidences
            accuNevents = accuNevents +...
                length(energyScanArchive{ii}.incidences{jj});
        end
        if accuNevents > 300000 || ii == nTrials
            chunkIdx    =   chunkIdx + 1;
            archiveChunk    =   energyScanArchive(startingIdx:ii);
            save(strcat(basepath,filename,'_',num2str(chunkIdx)),'archiveChunk')
            startingIdx = ii+1;
            accuNevents = 0;
        end
    end
    result = chunkIdx;
end