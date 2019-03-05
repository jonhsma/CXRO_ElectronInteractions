function [summary,delta]  = trajLengthAna(trajectory)
    % Compute the difference matrices
    r_ini = [trajectory.xyz_init];
    r_final = [trajectory.xyz_final];

    delta.M_ini = sqrt((r_ini(1,:)'-r_ini(1,:)).^2+...
        (r_ini(2,:)'-r_ini(2,:)).^2+...
        (r_ini(3,:)'-r_ini(3,:)).^2);

    delta.M_final = sqrt((r_final(1,:)'-r_final(1,:)).^2+...
        (r_final(2,:)'-r_final(2,:)).^2+...
        (r_final(3,:)'-r_final(3,:)).^2);
    
    % The results arrays. Thanks to the fact that two dimensional vector
    % indicing is not working
    deltaStepIni  = delta.M_final;
    deltaStepIni(:)   =   NaN;
    deltaStepFinal    =   deltaStepIni;
    
    v = ones([1 size(deltaStepIni,2)]);
    v(:) = NaN;
    summary.ini.mean    =   v;
    summary.ini.ms     =   v;
    summary.ini.std     =   v;
    summary.ini.sstd    =   v;
    summary.ini.N       =   v;
    
    summary.final = summary.ini;
    
    for jj = 1:length(trajectory)
        % The dummy vector
        v = ones([1 length(trajectory)]);
        v(:) = NaN;
        
        % The initial coordinates
        % Feeding the diagonals into the dummy
        u = diag(delta.M_ini ,jj);
        v(1:length(u)) = u;
        % Feeding the dummy into the output, which involves no possibility
        % of resizing the output
        % Initial coordinate difference sorted by the difference is step
        % numer
        deltaStepIni(jj,:)    =   v;
        % Statistics per speration
        summary.ini.mean(jj)    =   mean(u);
        summary.ini.ms(jj)      =   mean(u.^2); % This is actually stupid but since it doesn't slow things down too much
        summary.ini.std(jj)     =   std(u);
        summary.ini.sstd(jj)    =   std(u.^2);
        summary.ini.N(jj)       =   length(u);
        
        % The final coordinates
        u = diag(delta.M_final ,jj);
        v(1:length(u)) = u;
        % Feeding the dummy into the output, which involves no possibility
        % of resizing the output
        deltaStepFinal(jj,:)    =   v;
        % Statistics per speration
        summary.final.mean(jj)    =   mean(u);
        summary.final.ms(jj)     =   mean(u.^2);
        summary.final.std(jj)     =   std(u);
        summary.final.sstd(jj)    =   std(u.^2);
        summary.final.N(jj)       =   length(u);
        
    end
    summary.ini.deltaStep = deltaStepIni;
    summary.final.deltaStep = deltaStepFinal;
        
        
end