%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     pagidx      index of the pag that is the closest to the ref. pos.
%%%     npags       number of pags nearby (within the reaction radius)
%%%     polymidx    indices of ionizable polymers nearby
%%%     npolyms     number of ionizable polymers nearby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pagidx,npags,polymidx,npolyms]=pag_v_polym(xyz,pag_grid,pagimg,polym_img,rcnrad)

xgrid=pag_grid.x;
ygrid=pag_grid.y;
zgrid=pag_grid.z;

xref=xyz(1);
yref=xyz(2);
zref=xyz(3);

%% masking out the area covered by the reaction radius
tmpmask=zeros(size(xgrid));

% susmq=(xgrid-xref).^2+(ygrid-yref).^2;
% tmp=(susmq<=rcnrad^2).*sqrt(rcnrad^2-susmq);
% tmpmask(abs(zgrid-zref)<=tmp & tmp~=0)=1;

idx=find((xgrid-xref).^2+(ygrid-yref).^2+(zgrid-zref).^2<=rcnrad^2);
tmpmask(idx)=1;

if sum(sum(sum(tmpmask==1)))~=0
    dbg=1;
end

%% Identifying the usable pags and locate the closest one
pag_removable=pagimg.*tmpmask;
polym_ionizable=polym_img.*tmpmask;

pagidx=find(pag_removable~=0);
npags=sum(pag_removable(pagidx));
diff=[];
for i = 1:length(pagidx)
%     [xidx_tmp,yidx_tmp,zidx_tmp]=ind2sub(size(xgrid),pagidx(i));
    diff(i)=sqrt((xref-xgrid(pagidx(i))).^2+(yref-ygrid(pagidx(i))).^2+(zref-zgrid(pagidx(i))).^2);
end
if ~isempty(diff)
    tmpidx=find(diff==min(diff));
%     pagidx=pagidx(tmpidx(1));
    pagidx=pagidx(tmpidx);
end
%% Identifying the ionizable polymers
polymidx=find(polym_ionizable~=0);
npolyms=sum(polym_ionizable(polymidx));
end

%%%%%%%% defaults if no PAG is found within vicinity
% yesno=0;
% navail=sum(pag_removable(:));
% pagidx=idx;
%     
% %%%% if PAGs found within vicinity, update the return variables
% if sum(pag_removable(:))~=0
%     yesno=1;
%     navail=sum(pag_removable(:));
%     pagidx=idx;
% end
