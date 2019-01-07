%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function saves the database in chunks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = saveDatabase(simDatabase,basepath,filename)
    nTrials = length(simDatabase);
    accuNevents = 0;
    startingIdx = 1;
    chunkIdx    =   0;
    
    for ii = 1 : nTrials
        accuNevents = accuNevents + length(simDatabase(ii).eventsArray);
        if accuNevents > 300000 || ii == nTrials
            chunkIdx    =   chunkIdx + 1;
            databaseChunk    =   simDatabase(startingIdx:ii);
            save(strcat(basepath,filename,'_',num2str(chunkIdx)),'databaseChunk')
            startingIdx = ii+1;
            accuNevents = 0;
        end
    end
    result = chunkIdx;
end