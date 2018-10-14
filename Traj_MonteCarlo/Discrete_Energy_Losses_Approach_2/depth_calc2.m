function [meandepth,meandepth2,depthmat]=depth_calc2(deprimg,zvec,elecimg,thr)

for i = 1:size(deprimg,2)
    for j = size(deprimg,3):-1:1
        deprvec=deprimg(:,i,j);
        deprvec(deprvec>=thr)=1;
        deprvec(deprvec<thr)=0;
        idx=find(deprvec==1);
        if i==30 & j==30
            dbg=1;
        end
        if ~isempty(idx)
%             if idx(1)~=1
%                 depthmat(j,i)=0;
%             else
                depthmatB(j,i)=zvec(idx(end));
                depthmat(j,i)=sum(zvec.*deprvec)/length(deprvec(deprvec~=0));
%                 depthmat(j,i)=(length(idx)).*mean(diff(zvec));
%                 if depthmat(j,i)==4
%                     dbg=1;
%                 end
%             end
        else
            depthmat(j,i)=0;
            depthmatB(j,i)=0;
        end
    end
end

depthmat=depthmat(end:-1:1,:); % flip bottom-to-top to stay consistent with electron image co-ordinate system
[xidx,yidx]=find(elecimg>0);

depthmat_2=depthmat(min(xidx):max(xidx),min(yidx):max(yidx));
depthmat_2B=depthmatB(min(xidx):max(xidx),min(yidx):max(yidx));

% depthmat_2=depthmat;

idx=find(depthmat_2>0);
if ~isempty(idx)
%     meandepth=mean(depthmat_2(:)); % depthmat_2 is already the cropped one
    meandepth=mean(depthmat_2(idx)); % take only the non-zero depths. WORKED FOR 40 eV
    meandepth2=mean(depthmat_2B(idx)); 
%     meandepth=mean(depthmat(depthmat>0));
%     meandepth=mean(depthmat(:)); % shouldn't do this
else
    meandepth=0;
    meandepth2=0;
end