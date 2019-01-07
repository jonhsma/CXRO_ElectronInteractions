% To locate the least efficient line
tic
testDB = reorderArchive(energyScanArchive(1:2));
toc
% 12 seconds for 2
%% The inconsistency
% To investigate if the descrapency between counting energy discontinuity
% and summing nSE
size([testDB.eventsArray([testDB.eventsArray.nSE]>0)])
size([testDB.eventsArray([testDB.eventsArray.Ese]<2 &...
    [testDB.eventsArray.Ese]>0)])
%% Trying to figure out which discontinuity went wrong
flag_SEpending = 0;
nDiscon = 0;
for ii = 2:length(testDB.eventsArray)
    if abs(testDB.eventsArray(ii).Ein - testDB.eventsArray(ii-1).Eout) > 0 ||...
            testDB.eventsArray(ii).Ein == 85 % In the case of escaped primaries their life is one step long
        nDiscon = nDiscon + 1;
    end
end
disp(nDiscon);
%% Investigate the use of parfor to see if the tast is CPU intensive
tic
parfor ii = 1:22
    tesDBPar{ii}    =   reorderArchive(energyScanArchive(ii))
end
toc
% 24 seconds for 22 trials so it is CPU intensive
%% Storing things in an array. Is that a better idea?
tic
parfor ii = 1:22
    tesDBPar2(ii)    =   reorderArchive(energyScanArchive(ii))
end
toc
%% Cell and strubture arrays take up the same amount of space. Use array for easier access
tic
parfor ii = 1:88
    finalDB(ii)    =   reorderArchive(energyScanArchive(ii))
end
toc
%% Apparently saving is an ordeal. Try to see how I can export the database
fdb1 = finalDB(1:44);
save(strcat(outputBasePath,'finalFirstHalf'),'fdb1')
%% The last attempt worked but how far can I go?
fdb1 = finalDB(1:44);
save(strcat(outputBasePath,'finalFirstHalf'),'fdb1')

%% Write a save database function
saveDatabase(finalDB,outputBasePath,'DB_85eV')

%% And one to load the database
reconDB = loadDatabase(outputBasePath,'DB_85eV');

