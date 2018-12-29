function mfp = genMFP_Vibr(scattdataVibr,incidentE,moleculeDensity)
    E   =   scattdataVibr.ics(:,1);
    if incidentE > max(E)
        mfp = inf;
        return
    end
    %% Find the closest energy entry in the crosssection table
    idx1    =   find(E<incidentE);
    [~,idxtmp]=min(abs(E(idx1)-incidentE));
    idx1=idx1(idxtmp);

    idx2=find(E==incidentE);

    idx3=find(E>incidentE);
    [~,idxtmp]=min(abs(E(idx3)-incidentE));
    idx3=idx3(idxtmp);

    %% generate the ICS vector at desired energy
    %icsvec=[];
    if ~isempty(idx2) 
        Eochoose=idx2;
    else 
        Eochoose=[idx1 idx3];
    end
    Etmp=scattdataVibr.ics(:,1);
    ics2=scattdataVibr.ics(:,2:end);
    ics2(isnan(ics2))=0;
    ics2=sum(ics2,2);
    if length(Eochoose)==1
        icsval=ics2(Eochoose);
    else    
        icsval=lineInterpol([Etmp(Eochoose(1)) Etmp(Eochoose(2))],...
            [ics2(Eochoose(1)) ics2(Eochoose(2))],...
            incidentE);
    end
    % unit conversion from Angstrom^2 to cm^2
    icsTot=icsval*scattdataVibr.ics_mult;
    invImfp_invCM      =   icsTot*moleculeDensity;
    mfp = 1/invImfp_invCM*1e7;
end