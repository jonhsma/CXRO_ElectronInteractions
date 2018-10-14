function [defint,indefint]=integrate(x,y)

ysum=cumsum(y);
indefint=ysum(2:end).*(x(2:end)-x(1:end-1));

indefint=[];
for i = 2:length(y)
    indefint(i)=trapz(x(1:i),y(1:i));
end
indefint(1)=indefint(2);

defint=trapz(x,y);