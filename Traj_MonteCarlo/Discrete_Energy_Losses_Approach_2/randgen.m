function randval=randgen(x,CDF,N)
% fprintf('In randgen now\n');
% tic
randval=zeros(1,N);
for i = 1:N
    u=rand;
    idx=find(CDF==u);
    if ~isempty(idx)
        randval(i)=x(idx);
    else
        idx1=find(CDF<u);idx1=idx1(end);
        idx2=find(CDF>u);idx2=idx2(1);
        p=polyfit(x(idx1:idx2),CDF(idx1:idx2),1);
        randval(i)=(u-p(2))/p(1);
    end
end
% time=toc;
% fprintf('randgen: N = %d; took %.2f s = %.2f min = %.2f hours\n',N,time,time/60,time/3600);