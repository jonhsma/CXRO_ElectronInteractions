function simdata=gensimdata(dose,rt,rtlost_std)

simdata.dose=dose;
simdata.rt=rt;
simdata.rt_eb=rtlost_std;
% simdata.slope=polyfit(log10(simdata.dose),simdata.rt,1);
