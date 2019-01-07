%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function loads the database and stiches the chunks together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stichedDatabase = loadDatabase(basepath,filename)
    fileCounter = 1;
    stichedDatabase = [];
    while(exist(strcat(basepath,filename,'_',num2str(fileCounter),'.mat'),'file'))
        currChunk = load(strcat(basepath,filename,'_',num2str(fileCounter)));
        stichedDatabase =  [stichedDatabase currChunk.databaseChunk];
        fileCounter = fileCounter +1;
    end
end