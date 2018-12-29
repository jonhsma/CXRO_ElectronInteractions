function mfp = genMFP_OptData(eventInc,scattdataOptical,eIni)
    controlparm.onlyimfp=1;
    if isfield(scattdataOptical,'inel_dcsdata')
        Elossrand_opt=genrandEloss_OptData_JHM(scattdataOptical,eIni,scattdataOptical.inel_dcsdata,controlparm);
    else
        Elossrand_opt=genrandEloss_OptData_JHM(scattdataOptical,eIni,controlparm);
    end
    if eIni>eventInc.lowEthr
        mfp=Elossrand_opt.imfp; % imfp in nm; the above lines only calculate imfp due to the controlparms.onlyimfp line abvove.
        if isnan(mfp)
            mfp=Inf;
        end
    else
        mfp=eventInc.lowEimfp;
    end
end