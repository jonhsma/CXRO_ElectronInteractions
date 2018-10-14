function dydx=deriv(x,y)
nhood=1; % dydx(n)=polyfit(x(n-1:n+1),y(n-1:n+1),1);
algo=2;
if length(x)==2
    dydx=(y(2)-y(1))/(x(2)-x(1));
    dydx=[dydx(1) dydx];
else
    dydx=[];
    if algo==1 % center-difference
        dydx=(y(3:end)-y(1:end-2))./(x(3:end)-x(1:end-2));
        dydx=[dydx(1) dydx dydx(end)];
    end
    
    if algo==2 % linear fit out to neighborhood=nhood around center
        for idx = 1:length(y)
            leftidx=max([idx-nhood 1]);
            rightidx=min([idx+nhood length(y)]);
            p=polyfit(x(leftidx:rightidx),y(leftidx:rightidx),1);

            dydx(idx)=p(1);
        end
    end
    
    if algo==3 % simple look-forward difference
        dydx=(y(2:end)-y(1:end-1))./(x(2:end)-x(1:end-1));
        dydx=[dydx(1) dydx];
    end
end
dbg=1;