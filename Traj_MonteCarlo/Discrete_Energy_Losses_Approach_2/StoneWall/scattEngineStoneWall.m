%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is the stone-wall scattering engine
%%% It's an adaptation of the narasimhan approach. When an electron falls
%%% below a certain energ, it excites an acid and it's gone.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outcome     =   scattEngineStoneWall(eIncident,scattData)    
    if eIncident < scattData.stoneWall.CUTOFF
        outcome.imfp    =   scattData.stoneWall.IMFP;
        outcome.eLoss   =   eIncident;
        outcome.rxnR    =   scattData.stoneWall.ACID_REACTION_RADIUS;
    else
        outcome.imfp    =   Inf;
        outcome.eLoss   =   0;
        outcome.rxnR    =   0;
    end
    
    outcome.theta   =   0;
    outcome.phi     =   0;
end