function mfp = genMFP_OptData(eventInc,scattdataOptical,eIni)
    controlparm.onlyimfp=1;
    if isfield(scattdataOptical,'inel_dcsdata')
        Elossrand_opt=genrandEloss_OptData_JHM(scattdataOptical,eIni,scattdataOptical.inel_dcsdata,controlparm);
    else
        Elossrand_opt=genrandEloss_OptData_JHM(scattdataOptical,eIni,controlparm);
    end
    % This cut-off exist because there aren't any data below that
    % The low energy threshold option LOW_ENERGY_BEHAVIOUR_BOUNDARY is not
    % used becase it also determines when the vibrational mechanism kicks
    % in, thus always creating a singularity.
    if eIni> 20
        mfp=Elossrand_opt.imfp; % imfp in nm; the above lines only calculate imfp due to the controlparms.onlyimfp line abvove.
        if isnan(mfp)
            mfp=Inf;
        end
    else
        mfp=eventInc.lowEimfp;
    end
end