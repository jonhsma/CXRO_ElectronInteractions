function depth=depth_calc(zvec,depr,thr)
algo=1;

if algo==1
    idx=find(depr>thr);

    if ~isempty(idx)
        idxmin=min(idx);
        idxmax=max(idx);

        %%% the left edge:
        if idxmin>1
            p=polyfit(zvec(idxmin-1:idxmin),depr(idxmin-1:idxmin),1);
            leftedge=(thr-p(2))/p(1);
        else
            leftedge=zvec(1);
        end

        %%% the right edge:
        if idxmax<length(zvec)
            p=polyfit(zvec(idxmax:idxmax+1),depr(idxmax:idxmax+1),1);
            rightedge=(thr-p(2))/p(1);
        else
            rightedge=zvec(end);
        end

        depth=rightedge-leftedge;
    else
        depth=0;
    end
end
