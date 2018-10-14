function outdata=getNthruz(img)
%%% usage: outdata=getNthruz(img)
%%% img: (INPUT) the 3-D image containing data
%%% outdata: (OUTPUT) the output data structure

N=zeros(size(img,1),1);
for i = 1:size(img,1)
    tmpmat=img(i,:,:);
    N(i,1)=sum(tmpmat(:));
    len=length(find(tmpmat(:)>0));
    N_per_px(i,1)=N(i,1)/len;
end

outdata.N=N;
outdata.N_per_px=N_per_px;