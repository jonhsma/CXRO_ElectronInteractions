function [xout,yout,yout_eb]=getMean(x,y,y_eb)

[xunique,idx]=unique(x);

xout=[];
yout=[];
for i = 1:length(xunique)
    idx=find(x==xunique(i));
    xout(i)=xunique(i);
    yout(i)=mean(y(idx));
    yout_eb(i)=mean(y_eb(idx));
end

