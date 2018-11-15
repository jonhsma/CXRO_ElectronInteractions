%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function displays the hierachy of the electrons of a certain
%%% incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function graphHandle = displayTrajHierachy(incidence)
    trajCollections = splitIncidence(incidence);
    primaryTrc = displayTrajectory(trajCollections{1},[0.25 1 1]);
    for jj = 2:length(trajCollections)
        secondaryTrc = displayTrajectory(trajCollections{jj},[0.75 0 0]);
    end
    if length(trajCollections) >=2
        legend([primaryTrc secondaryTrc],'Primary','Secondary')
    end
    graphHandle = gcf;
end