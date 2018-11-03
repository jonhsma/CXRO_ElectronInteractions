%%% This function standardize the acid removal process. When changes are
%%% made to acid activation, they need to be made once and right here.
function [posPAG,posPAG_removed,acid_act_xyz_idx,acid_act_e_xyz] =...
    acidActivation(num2Remove,posPAG,posPAG_removed,...
    pagidx,eventPosition,...
    acid_act_xyz_idx,acid_act_e_xyz)
    
    %%% Determine which pag to remove 
    remove_idx      =   randi([1 length(pagidx)],num2Remove);
    posPAG_removed  =   [posPAG_removed,...
        posPAG(:,pagidx(remove_idx))];
    posPAG(:,pagidx(remove_idx)) = NaN;

    %%%Register the acid generation events
    acid_act_xyz_idx    =   [acid_act_xyz_idx pagidx(remove_idx)];
    acid_act_e_xyz      =   [acid_act_e_xyz;...
        ones([length(remove_idx) 1])*eventPosition];
end