function outparms=RcnDiffn(physparms,simparms,data)
% if ~exist('poissrnd')
%     poissrnd=@randpois; % due to weird error in this version of 2013a somehow doesn't have randpois, so re-direct to PPN's implementation
% end
%%% extract the physical parameters
base_load=physparms.base_load;
prot_load=physparms.prot_load;
acid_difflen=physparms.acid_difflen;
base_difflen=physparms.base_difflen;
kD=physparms.kD;
kQ=physparms.kQ;

%%% extract the simulation parameters:
time=simparms.sim_time;
dt=simparms.dt_PEB;
rng_shuffle=simparms.rng_shuffle;

acidimg=data.acidimg;
% acidimg=10.*acidimg; % artificial dose increase

if rng_shuffle==1
    rng('shuffle');
else
    rng(102); % just a random integer here if we want deterministic randomness for concentrations
end
baseimg=base_load.*ones(size(acidimg)).*prod(data.univ.px_nm);
rng('shuffle');
baseimg=poissrnd(baseimg);
rng('default');

protimg=prot_load.*ones(size(acidimg)).*prod(data.univ.px_nm);
rng('shuffle');
protimg=poissrnd(protimg);
rng('default');

%%%% Add zeros in all of the image matrices
acidimg=[zeros(simparms.zadd_npx,size(acidimg,2),size(acidimg,3));acidimg];
baseimg=[zeros(simparms.zadd_npx,size(baseimg,2),size(baseimg,3));baseimg];
protimg=[zeros(simparms.zadd_npx,size(protimg,2),size(protimg,3));protimg];

protimg_pre=protimg;
rng('shuffle');

%%%% define the 3D PSF:
xgrid=data.univ.grid.x;
ygrid=data.univ.grid.y;
zgrid=data.univ.grid.z;
Acidsig=acid_difflen/sqrt(length(time));
Basesig=base_difflen/sqrt(length(time));

psf_radius=sqrt(xgrid.^2+ygrid.^2+zgrid.^2);
xyz_idx=find(psf_radius==min(min(min(psf_radius))));
[zo,xo,yo]=ind2sub(size(psf_radius),xyz_idx);
xo=xo(1);yo=yo(1);zo=zo(1); % just pick the first place where min occurs

% xo=-0.5;yo=-0.5;zo=-0.5;
xo=0;yo=0;zo=0; % may need to uncomment above line if you see the deprot. image walking or something
xo=xgrid(23,23,23);yo=ygrid(23,23,23);zo=zgrid(23,23,23);

nz=size(xgrid,1);nx=size(xgrid,2);ny=size(xgrid,3);

if mod(nz,2)==0 % is it even?
    xo=mean(mean(mean(xgrid(nz/2:nz/2+2,nz/2:nz/2+2,nz/2:nz/2+2))));
    yo=mean(mean(mean(ygrid(nz/2:nz/2+2,nz/2:nz/2+2,nz/2:nz/2+2))));
    zo=mean(mean(mean(zgrid(nz/2:nz/2+2,nz/2:nz/2+2,nz/2:nz/2+2))));
else
    xo=xgrid(floor(nz/2)+1,floor(nz/2)+1,floor(nz/2)+1);
    yo=ygrid(floor(nz/2)+1,floor(nz/2)+1,floor(nz/2)+1);
    zo=zgrid(floor(nz/2)+1,floor(nz/2)+1,floor(nz/2)+1);
end
psf_radius=sqrt((xgrid-xo).^2+(ygrid-yo).^2+(zgrid-zo).^2);

acidpsf=1/(Acidsig^3*(2*pi)^(3/2)).*exp(-psf_radius.^2./(2*Acidsig^2));
basepsf=1/(Basesig^3*(2*pi)^(3/2)).*exp(-psf_radius.^2./(2*Basesig^2));
acidpsf=acidpsf./sum(acidpsf(:)); % normalize as sum
basepsf=basepsf./sum(basepsf(:)); % normalize as sum


%%%% develop a heuristic to crop the PSF size so the convolution runs faster
numsig=3;  % number of sigmas to keep [cropping the rest]
% numsig_npx=ceil(numsig/prod(data.univ.px_nm)); % Looks WRONG

numsig_npx=ceil(numsig*Acidsig/(data.univ.px_nm(1))); % Be wary about the px_nm being a vector, for the future.
xsel_idx=max([1,floor(size(acidpsf,2)/2-numsig_npx)]):min([floor(size(acidpsf,2)),floor(size(acidpsf,2)/2+numsig_npx)]);
ysel_idx=max([1,floor(size(acidpsf,3)/2-numsig_npx)]):min([floor(size(acidpsf,3)),floor(size(acidpsf,2)/2+numsig_npx)]);
zsel_idx=max([1,floor(size(acidpsf,1)/2-numsig_npx)]):min([floor(size(acidpsf,1)),floor(size(acidpsf,2)/2+numsig_npx)]);
% acidpsf=acidpsf./sum(acidpsf(:)); % normalize as sum
acidpsf=acidpsf(zsel_idx+1,xsel_idx+1,ysel_idx+1);
% basepsf=basepsf./sum(basepsf(:)); % normalize as sum
basepsf=basepsf(zsel_idx+1,xsel_idx+1,ysel_idx+1);

% acidpsf=acidpsf(1:end-1,1:end-1,1:end-1); % may prevent the deprot. image walking issue
% basepsf=acidpsf(1:end-1,1:end-1,1:end-1); % may prevent the deprot. image walking issue

numsig_npx=ceil(numsig*Basesig/(data.univ.px_nm(1))); % Be wary about the px_nm being a vector, for the future.
xsel_idx=max([1,floor(xo-numsig_npx)]):min([size(basepsf,2),ceil(xo+numsig_npx)]);
ysel_idx=max([1,floor(yo-numsig_npx)]):min([size(basepsf,3),ceil(yo+numsig_npx)]);
zsel_idx=max([1,floor(zo-numsig_npx)]):min([size(basepsf,1),ceil(zo+numsig_npx)]);
% acidpsf=acidpsf./sum(acidpsf(:)); % normalize as sum
% basepsf=basepsf(zsel_idx,xsel_idx,ysel_idx);

%%%% convert rcn rates to pixels/s
kD_px=kD/prod(data.univ.px_nm);
kQ_px=kQ/prod(data.univ.px_nm);

dbg=1;

for j = 1:length(time)
    fprintf('...RcnDiffn: PEB Time %d of %d\n',j,length(time));
    tstart_PEBStep=tic;
    
    %%%% Diffusion of acids/bases
    acidimg=convn(acidimg,acidpsf,'same');
    if base_difflen~=0
        baseimg=convn(baseimg,basepsf,'same');
    end
    
    %%%% Acid/Base and deprotection reactions:
%     conc_ratio=baseimg./(baseimg+protimg); % relative concentrations
%     conc_ratio(baseimg==0)=0;
%     conc_ratio2=protimg./(baseimg+protimg);
%     conc_ratio2(protimg==0)=0;
    
    conc_ratio=1;
    conc_ratio2=1; % the above sometimes make quenching much sloer than deprotection when amount of base is tiny, not correct. this is better
    
    %%%%%%% update the matrices, subtracting the appropriate values:
    %%% Acid/Base reaction:
    Nrcn_AB=dt*kQ_px.*acidimg.*baseimg.*conc_ratio;
    acidimg=acidimg-Nrcn_AB;acidimg(acidimg<0)=0;
    baseimg=baseimg-Nrcn_AB;baseimg(baseimg<0)=0;

    %%% Deprotection reaction:
    Ndeprot=dt*kD_px.*acidimg.*protimg.*conc_ratio2;
%     Ndeprot=dt*kD_px.*acidimg.*protimg;
    protimg=protimg-Ndeprot;protimg(protimg<0)=0;
    
    tend_PEBStep=toc(tstart_PEBStep);
    fprintf('......Took %.4f s\n',tend_PEBStep);
end
deprimg=(protimg_pre-protimg)./protimg_pre;
deprimg(isnan(deprimg))=0;
deprimg(deprimg==Inf)=0;
deprimg(deprimg==-Inf)=0;
% deprimg(protimg_pre==0)=1; % if no prot. group, it will dissolve [that was stupid!]

%%%% assign the output variables:
outparms.deprimg=deprimg;
outparms.protimg_pre=protimg_pre;
outparms.acidpsf=acidpsf;
outparms.basepsf=basepsf;