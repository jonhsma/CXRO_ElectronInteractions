%%% This script is a adaptaion to modularization of scattering mechanisms.
%%% It uses the engine and only spits out the MFP. I don't simplify it
%%% because the engine is simple enough that it shouldn't have been slowing
%%% the system down
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mfp = genMFP_StoneWall(eIncident,scattData)
    results = scattEngineStoneWall(eIncident,scattData);
    mfp = results.imfp;
end