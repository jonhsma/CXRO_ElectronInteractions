%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function initializes the optical scatteing data by integrating the
%%% integral cross-section ahead of time instead of on-demand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function optdata = optDataScatt_Init(optdata)
    dcsArray = optdata.inel_dcsdata;
    
    thetaScale      =   dcsArray.thetamat; % just a single vector in the most recent one
    dsigdOmega      =   dcsArray.dsigdOmega;
    
    
    isc = zeros(size(thetaScale));
    for jj = 1:size(dsigdOmega,1)
        for i = 2:size(dsigdOmega,2)
            isc(jj,i)=trapz(thetaScale(jj,1:i),dsigdOmega(jj,1:i).*sin(thetaScale(jj,1:i)));
        end
        isc(jj,:)     =   isc(jj,:)/isc(jj,end);
    end
    
    
    
    optdata.inel_dcsdata.isc = isc;
end

