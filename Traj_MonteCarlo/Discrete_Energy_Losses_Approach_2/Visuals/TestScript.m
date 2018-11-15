%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is a test script on the trajectory visuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nov 14 The individual trajectory display
testTraj = scanArchive{1}{46}.incidences{2};
displayTrajectory(testTraj,[1 1 1])
hold off
displayTrajectory(testTraj,'energy',[],[],30)

%% Actually there is a secondary. Test the secondary display
trajCollections = splitIncidence(testTraj);
hold off
displayTrajectory(trajCollections{1},[0.25 1 1])
displayTrajectory(trajCollections{2},[0.75 0 0])
%%% But it travelled too short of a distance

%% Get some 80 eV data to see secondaries
testTraj = scanArchive{1,52}.incidences{2};
figure(98000)
hold off
displayTrajectory(testTraj,'energy',[],[],80,jet(100))
colormap(jet)
caxis([0 80])
colorbar

figure(98001)
trajCollections = splitIncidence(testTraj);
hold off
displayTrajectory(trajCollections{1},[0.25 1 1])
for jj = 2:length(trajCollections)
    displayTrajectory(trajCollections{jj},[0.75 0 0])
end