function y = tabulateSummary(resultObj)

    acidSpread          =   zeros([size(resultObj,2) 5]);
    excitationSpread    =   zeros([size(resultObj,2) 5]);
    
    for ii = 1:size(resultObj,2)
        % First three columns are the standard diviations
        % The forth is the root of sum of the 3 sds
        % The fifth is the mean of r
        acidSpread(ii,1:3)          =   resultObj{ii}.sig_acid;
        excitationSpread(ii,1:3)    =   resultObj{ii}.sig_acid_fine;
        
        acidSpread(ii,4)    =   sqrt(sum(acidSpread(ii,1:3).^2));
        acidSpread(ii,5)    =   mean(resultObj{ii}.rabs);
        
        excitationSpread(ii,4)  =   sqrt(sum(excitationSpread(ii,1:3).^2));
        excitationSpread(ii,5)  =   mean(resultObj{ii}.rabs_fine);
    end
    disp('Acid positions')
    fprintf('%12s%12s%12s%12s%12s\n',...
        'sigma_x','sigma_y','sigma_z','(<r^2>)^0.5','<r>')
    for ii = 1:size(resultObj,2)
        fprintf('%12.3f%12.3f%12.3f%12.3f%12.3f\n',...
            acidSpread(ii,1),acidSpread(ii,2),acidSpread(ii,3),...
            acidSpread(ii,4),acidSpread(ii,5))
    end
    disp('Excitation positions')
    fprintf('%12s%12s%12s%12s%12s\n',...
        'sigma_x','sigma_y','sigma_z','(<r^2>)^0.5','<r>')
    for ii = 1:size(resultObj,2)
        fprintf('%12.3f%12.3f%12.3f%12.3f%12.3f\n',...
            excitationSpread(ii,1),excitationSpread(ii,2),excitationSpread(ii,3),...
            excitationSpread(ii,4),excitationSpread(ii,5))
    end
    
    y.acid          =   acidSpread;
    y.excitation    =   excitationSpread;

end