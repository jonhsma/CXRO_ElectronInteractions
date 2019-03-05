% This function returns the cosine between every pair in a trajectory
function [summary,dpM] =trajCorr(trajectory)
    % The correlation between all steps
    direction = [trajectory.xyz_delta]./sum([trajectory.xyz_delta].^2,1).^0.5;
    dpM = direction'*direction;
    
    % claim a block of memory
    dpM_byDist      =   ones(size(dpM));
    dpM_byDist(:,:)   =   NaN;
    
    % The results arrays. Thanks to the fact that two dimensional vector
    % indicing is not working
    v = ones([1 size(dpM,2)]);
    v(:) = NaN;
    summary.mean    =   v;
    summary.std     =   v;
    
    for jj = 1:length(trajectory)
        % The dummy vector
        v = ones([1 size(dpM,2)]);
        v(:) = NaN;
        % Feeding the diagonals into the dummy
        u = diag(dpM,jj);
        v(1:length(u)) = u;
        % Feeding the dummy into the output, which involves no possibility
        % of resizing the output
        dpM_byDist(jj,:)    =   v;
        % Statistics per speration
        summary.mean(jj)    =   mean(u);
        summary.std(jj)     =   std(u);
    end
    
    rmNaN = ~isnan(dpM_byDist);    
    summary.dp_byDist = dpM_byDist;    
    summary.N = sum(rmNaN,2);
end