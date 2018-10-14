function [yesno,navail,pagidx,npixels_vol]=pag_avail(xyz,pag_grid,pagimg,rcnrad)

xgrid=pag_grid.x;
ygrid=pag_grid.y;
zgrid=pag_grid.z;

xref=xyz(1);
yref=xyz(2);
zref=xyz(3);

tmpmask=zeros(size(xgrid));

% susmq=(xgrid-xref).^2+(ygrid-yref).^2;
% tmp=(susmq<=rcnrad^2).*sqrt(rcnrad^2-susmq);
% tmpmask(abs(zgrid-zref)<=tmp & tmp~=0)=1;

idx=find((xgrid-xref).^2+(ygrid-yref).^2+(zgrid-zref).^2<=rcnrad^2);
tmpmask(idx)=1;
npixels_vol=length(find(tmpmask==1));

if sum(sum(sum(tmpmask==1)))~=0
    dbg=1;
end

pag_removable=pagimg.*tmpmask;
idx=find(pag_removable~=0);

%%%%%%%% defaults if no PAG is found within vicinity
yesno=0;
navail=sum(pag_removable(:));
pagidx=idx;
    
%%%% if PAGs found within vicinity, update the return variables
if sum(pag_removable(:))~=0
    yesno=1;
    navail=sum(pag_removable(:));
    pagidx=idx;
end
