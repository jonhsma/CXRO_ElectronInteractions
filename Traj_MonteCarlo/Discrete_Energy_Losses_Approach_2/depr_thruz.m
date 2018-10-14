function deprdata=depr_thruz(deprdata_in,elecimg,acidimg,deprthr)
%%% usage: depr=depr_thruz(deprimg,elecimg,depr_thr)
%%% deprimg: (INPUT) the 3-D deprotection matrix
%%% elecimg: (INPUT) the electron 2-D image that shows the incident map
%%% depr_thr: (INPUT) the deprotection threshold
%%% outdata: (OUTPUT) the thru-z deprotection (mean and sigma)

deprimg=deprdata_in.deprimg;
protimg_pre=deprdata_in.protimg_pre;

[xidx,yidx]=find(elecimg>0);

xidx=[20 30]; % when you want to force the pixels where to look 
yidx=[20 30]; % when you want to force the pixels where to look 

for i = 1:size(deprimg,1)
    tmpmat=deprimg(i,:,:);
    tmpmat=reshape(tmpmat,size(tmpmat,2),size(tmpmat,3));
    
    protmat=protimg_pre(i,:,:);
    protmat=reshape(protmat,size(protmat,2),size(protmat,3));
    
    tmpmat2=tmpmat(min(xidx):max(xidx),min(yidx):max(yidx));
    tmpmat3=protmat(min(xidx):max(xidx),min(yidx):max(yidx));
    
    tmpmat4=tmpmat2;
    tmpmat4(tmpmat3==0)=1;
    
    deprmean3(i)=mean(tmpmat4(:));
    
    if any(tmpmat2(:)>0)
        deprmean2(i)=mean(tmpmat2(tmpmat2>0)); % everywhere >0, take mean; will avoid using this in the future
    else
        deprmean2(i)=0;
    end
    if isnan(deprmean2(i))
        dbg=1;
    end
    
    npixels(i)=length(find(tmpmat2(:)>0));
    
    if deprmean2(i)>deprthr
        tmpmat2(tmpmat3==0)=1; % where no protecting group, treat as soluble.
    end
    
    deprmean(i)=mean(tmpmat2(:));
    deprstd(i)=std(tmpmat2(:));
end

deprdata.mean=deprmean;
deprdata.mean2=deprmean2;
deprdata.mean3=deprmean3;
deprdata.npxgood=npixels;
deprdata.std=deprstd;