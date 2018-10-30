%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     pagidx      index of the pag that is the closest to the ref. pos.
%%%     npags       number of pags nearby (within the reaction radius)
%%%     polymidx    indices of ionizable polymers nearby
%%%     npolyms     number of ionizable polymers nearby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     xyz         [x y z]
%%%     posPAG(:,i) [...,[x;y;z],...] of the i-th PAG and the same goes for
%%%                 the polymers
function [pagIdx,nPags,polymIdx,nPolyms]=pag_polym_query(xyz,posPAG,posPolymer,rxnRad)

posRef = xyz';

%% Locating the PAGs and polymers within the reaction rdius
%%% distance of each pag from reference 
distPAG     =   sqrt(sum((posPAG - posRef).^2,1));
idxPAG      =   find(distPAG < rxnRad);
distPolymer =   sqrt(sum((posPolymer - posRef).^2,1));
idxPolymer  =   find(distPolymer < rxnRad);
%% Identifying the usable pags and locate the closest one(s)
nPags       =   length(idxPAG);
%%% retrive the distances of pags to the reference position for PAGs within
%%% the reaction radius
diff_sub    =   distPAG(idxPAG);

if ~isempty(diff_sub)
    %%% Returning the indicies (for quering posPAG array) of the relavant PAGs 
    pagIdx  =   idxPAG(diff_sub==min(diff_sub));
end
%% Identifying the ionizable polymers
%%% Total number of polymers inside the reaction radius
nPolyms     =   length(idxPolymer);
polymIdx    =   idxPolymer;
end