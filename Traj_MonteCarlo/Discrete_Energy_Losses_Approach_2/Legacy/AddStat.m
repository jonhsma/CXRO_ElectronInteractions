function outdata=AddStat(indata,adddata)
%%% usage: outdata=AddStat(indata)
%%% indata: (INPUT) Data structure where data is to be added (dose vs. rt)
%%% outdata: (OUTPUT) indata, but the new data set added

dose=indata.dose;
rt=indata.rt;

dose2=adddata.dose;
rt2=adddata.rt;

for i = 1:length(dose2)
    idx=find(dose==dose2(i));
    if ~isempty(idx)
        indata.rt(idx)=mean([indata.rt(idx) rt2(i)]);
    else
        indata.dose(end+1)=dose2(i);
        indata.rt(end+1)=rt2(i);
    end
end

outdata=indata;