function outparms=Count_thruz(img,elecimg)
depr_thr=0.5;
thr_min=0.01;

[xidx,yidx]=find(elecimg>0);
xsel_lims=[min(xidx) max(xidx)];
ysel_lims=[min(yidx) max(yidx)];

for zcount=1:size(img,1) % go through each slice in z-coordinate, get the number of acids
    slice=img(zcount,:,:);
    slice=reshape(slice,[size(slice,2) size(slice,3)]);
    slice(slice<thr_min)=0;
    
    nquanta(zcount)=sum(slice(:));
    
    if all(all(slice==0)==1)
        slice_mean(zcount)=0;
    else
%         slice_mean(zcount)=mean(slice(slice>0));
        slice_mean(zcount)=mean(mean(slice(xsel_lims(1):xsel_lims(2),ysel_lims(1):ysel_lims(2))));
    end

end

outparms.N=nquanta;
outparms.mean=slice_mean;