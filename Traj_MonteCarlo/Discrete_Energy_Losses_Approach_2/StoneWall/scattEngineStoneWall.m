%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is the stone-wall scattering engine
%%% It's an adaptation of the narasimhan approach. When an electron falls
%%% below a certain energ, it excites an acid and it's gone.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outcome     =   scattEngineStoneWall(eIncident,scattDataStoneWall)    
    if eIncident < scattDataStoneWall.CUTOFF
        outcome.imfp    =   scattDataStoneWall.IMFP;
        outcome.eLoss   =   eIncident+5; % + 5 to ensure termination, even if trajectory terminates at zero
        outcome.rxnR    =   scattDataStoneWall.ACID_REACTION_RADIUS;
    else
        outcome.imfp    =   Inf;
        outcome.eLoss   =   0;
        outcome.rxnR    =   0;
    end
    
    outcome.theta   =   0;
    outcome.phi     =   0;
end