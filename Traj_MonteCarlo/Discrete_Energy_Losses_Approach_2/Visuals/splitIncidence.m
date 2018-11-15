%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function splits an incidence with secondaries and give each
%%% electron its only trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trajectories = splitIncidence(events)
    
    currentTraj = 1;
    
    trajectories = cell(1);
    trajectories{currentTraj}   =   events(1);
    for ii = 2:length(events)
        if sum((events{ii}.xyz_init-events{ii-1}.xyz).^2) <=0.001
            trajectories{currentTraj} = [trajectories{currentTraj} events(ii)];
        else
            currentTraj     = currentTraj+1;
            trajectories    = [trajectories cell([1 1])];
            trajectories{currentTraj} = [trajectories{currentTraj} events(ii)];
        end
    end

end