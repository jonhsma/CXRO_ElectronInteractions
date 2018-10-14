%% LEEM simulations analysis: Load files
clear;
clc;
close all;

% [fname_tmp,pathname]=uigetfile(['LEEM_SimResults_2\*.mat'],'MultiSelect','on');
% [fname_tmp,pathname]=uigetfile(['LEEM_Sims_ThruTrial\Ein=91*.mat'],'MultiSelect','on');
% [fname_tmp,pathname]=uigetfile(['LEEM_Sims_ThruTrial\*.mat'],'MultiSelect','on');
[fname_tmp,pathname]=uigetfile(['LEEM_SmallestRcnRad\*.mat'],'MultiSelect','on');

if ~iscell(fname_tmp)
    fnames{1}=(fname_tmp);
else
    fnames=fname_tmp;
end

%% LEEM simulations analysis: Read the loaded files
clc;

simdata={};

elecimg_bin=[];
for i = 1:length(fnames)
%     tmp=load([pathname fnames{i}]);
    tmp=load([pathname fnames{i}],'acidimg','zval','Dose','univ','Esweep','E_count','nacids_thrutrial','nacids_unsat_thrutrial','nacids_unsat_total');
    simdata{i}.E=tmp.Esweep(tmp.E_count);
    simdata{i}.acidimg=tmp.acidimg;
    simdata{i}.univ=tmp.univ;
%     simdata{i}.elecimg_inc=tmp.elecimg_inc;
%     simdata{i}.Dose=tmp.Dose;
    simdata{i}.zval=tmp.zval;
    simdata{i}.nAcids=sum(tmp.acidimg(:));
    
    %%% the 2 lines below are for sat/unsat acids comparisons
    simdata{i}.nacids_thrutrial=tmp.nacids_thrutrial;
    simdata{i}.nacids_unsat_thrutrial=tmp.nacids_unsat_thrutrial;
    simdata{i}.nacids_unsat_total=tmp.nacids_unsat_total;
    %%% end of the 2 lines for sat/unsat acids comparisons
    
    simdata{i}.Dose=tmp.Dose;
    if isfield(tmp,'elecimg_bin')
        if isempty(elecimg_bin)
            elecimg_bin=tmp.elecimg_bin;
        else
            elecimg_bin=elecimg_bin+tmp.elecimg_bin;
        end
    end
    tmp=[];
end
elecimg_bin=elecimg_bin>0;

%% post-proc of LEEM analysis results: RCN/Diffn/Base Model [set up parameters]
clc;
% close all;

base_load=0.04; % PAG is 0.5/nm3;
acid_difflen=8; % nm, over a 90 s PEB time at 100 deg. C.
base_difflen=5; % nm, over a 90 s PEB time at 100 deg. C.
prot_load=1; % /nm3
kD=1; % Deprotection Rate (nm3/s)
kQ=2; % A/B Quenching Rate (nm3/s)
tPEB=90;
dt_PEB=1;
time=dt_PEB:dt_PEB:tPEB;

physparms.base_load=base_load;
physparms.acid_difflen=acid_difflen;
physparms.base_difflen=base_difflen;
physparms.prot_load=prot_load;
physparms.kD=kD;
physparms.kQ=kQ;

simparms.tPEB=tPEB;
simparms.dt_PEB=dt_PEB;
simparms.sim_time=time;
simparms.rng_shuffle=0; % 0: use deterministic randomness; 1: shuffle

% fig1=figure;
% hold on;

% fig2=figure;
% hold on;

plt_styles={'-b','-r','-k','-g','--b','--r','--k','--g'};
plt_styles=repmat(plt_styles,[1,100]);

zvec=[];
deprmean=[];
depth=[];
depth2=[];
dose=[];

%% post-proc of LEEM analysis results: RCN/Diffn/Base Model [Just get the rcn/diffn psf]
clc;
close all;

tmpdata=simdata{1}; % just a placeholder to get grids etc.
tmpdata.acidimg=zeros(tmpdata.univ.npx);
[nx,ny,nz]=size(tmpdata.acidimg);
xyz_tmp=round([nx,ny,nz]./2);
xo=xyz_tmp(1);yo=xyz_tmp(2);zo=xyz_tmp(3);
xo=25;yo=25;zo=25;
tmpdata.acidimg(xo,yo,zo)=1;

simparms.zadd_npx=0;

outparms=RcnDiffn(physparms,simparms,tmpdata);

clc;
close all;

deprimg=outparms.deprimg;
% deprimg=deprimg./max(deprimg(:));
depr_zx=sum(deprimg,3);
depr_zx=depr_zx./max(depr_zx(:));

depr_xy=reshape(sum(deprimg,1),size(deprimg,2),size(deprimg,3));

tmpidx=find(deprimg>=0.5);
[xidx,yidx,zidx]=ind2sub(size(deprimg),tmpidx);
radius=sqrt((xidx-xo).^2+(yidx-yo).^2+(zidx-zo).^2);

%%% versus radius
xo=25;yo=25;
dx=1;dy=1;
x=1:dx:size(deprimg,1);
[x,y]=meshgrid(x);
x=x-xo;y=y-yo;
r=sqrt(x.^2+y.^2);

rvec=r(:);
deprvec=depr_zx(:);

deprmean=[];    deprstd=[];

rvec_unique=unique(rvec);rvec_unique=rvec_unique';
for j = 1:length(rvec_unique)
    idx=find(rvec==rvec_unique(j));
    deprmean(1,j)=mean(deprvec(idx));
    deprstd(1,j)=std(deprvec(idx));
end

figure;hist(radius(:),50)

% figure;imagesc(sum(deprimg,3));colorbar;

figure; imagesc(depr_zx); title('zx'); colorbar;

figure; contour(depr_zx,[0.5]); title('zx contour'); colorbar;

% figure; imagesc(depr_xy./max(depr_xy(:))); title('xy'); colorbar;

% figure;
% plot(depr_zx(25,:),'-ob','linewidth',3.0,'MarkerSize',12);
% hold on;
% title('zx cross-section');

figure;
plot(rvec_unique,deprmean(1,:),'-b','linewidth',3.0);

%% analyze saved rcn/diffusion psf files
clear;
clc;
close all;

data1=load('RcnDiffn_PSFData\PSF3D_Bblur=5_Ablur=26_kD=1.5_kQ=3.mat');
data2=load('RcnDiffn_PSFData\PSF3D_Bblur=5_Ablur=12_kD=1_kQ=2.mat');

dr=1;
radmin=0;
radmax=60;
radvec=radmin:dr:radmax;

xo=25;yo=25;
dx=1;dy=1;
x=1:dx:size(data1.deprimg,1);
[x,y]=meshgrid(x);
x=x-xo;y=y-yo;
r=sqrt(x.^2+y.^2);

rvec=r(:);
deprvec_1=data1.depr_zx(:);
deprvec_2=data2.depr_zx(:);

deprmean=[];    deprstd=[];

rvec_unique=unique(rvec);rvec_unique=rvec_unique';
for j = 1:length(rvec_unique)
    idx=find(rvec==rvec_unique(j));
    deprmean(1,j)=mean(deprvec_1(idx));
    deprstd(1,j)=std(deprvec_1(idx));
end

for j = 1:length(rvec_unique)
    idx=find(rvec==rvec_unique(j));
    deprmean(2,j)=mean(deprvec_2(idx));
    deprstd(2,j)=std(deprvec_2(idx));
end

figure;
imagesc(data1.depr_zx);colorbar;title('Ablur=24nm');

figure;
imagesc(data2.depr_zx);colorbar;title('Ablur=12nm');

figure;
plot(data1.depr_zx(25,:),'b','linewidth',3.0);
hold on;
plot(data2.depr_zx(25,:),'k','linewidth',3.0);
set(gca,'fontsize',30,'linewidth',3.0);

figure;
plot(r(:),data1.depr_zx(:)./max(data1.depr_zx(:)),'-b')
hold on;
plot(r(:),data2.depr_zx(:)./max(data2.depr_zx(:)),'-r')

figure;
errorbar(rvec_unique,deprmean(1,:),deprstd(1,:),'-b','linewidth',2.0);
hold on;
errorbar(rvec_unique,deprmean(2,:),deprstd(2,:),'-r','linewidth',2.0);
ylim([0 1.2]);
% xlim([-0.1 15]);


figure;
plot(rvec_unique,deprmean(1,:)./deprmean(1,1),'-b','linewidth',3.0);
hold on;
plot(rvec_unique,deprmean(2,:)./deprmean(2,1),'-r','linewidth',3.0);
ylim([0 1.1]);
% xlim([-0.5 20]);
set(gca,'fontsize',30,'linewidth',3.0);

%% post-proc of LEEM analysis results: RCN/Diffn/Base Model [Run the actual sim]
% close all;

dose=[];
zvec2=[];

ntrials=1;
deprthr=0.5;

deprmean=[];deprmean2=[];deprmean3=[];
deprstd=[];deprstd2=[];deprstd3=[];
tstart=tic;
fig1=figure;
simparms.zadd_npx=0;
dose=zeros(1,length(simdata));
rcndiffndata={};
for i = 1:length(simdata)
% for i = 7
    deprmean_tmp=[];deprstd_tmp=[];
    deprmean2_tmp=[];deprstd2_tmp=[];
    deprmean3_tmp=[];deprstd3_tmp=[];
    for trial_count=1:ntrials
    %     tloop_start=tic;
        waitbar(((i-1)*ntrials+trial_count)/(length(simdata)*ntrials));
        fprintf('Sweep simdata %d of %d\n',i,length(simdata));
        dose(1,i)=simdata{i}.Dose;
%         simdata_tmp=simdata{i};     simdata_tmp.acidimg=simdata_tmp.acidimg(1:25,:,:);
        outparms=RcnDiffn(physparms,simparms,simdata{i});
%         outparms=RcnDiffn(physparms,simparms,simdata_tmp);

        rcndiffndata{i}=outparms;
        
        figure(fig1);
        imagesc(mean(outparms.deprimg(:,:,20:30),3));caxis([0 1]);colorbar;
        drawnow;
    %     NdeprData=Count_thruz(outparms.deprimg,simdata{i}.elecimg_inc);
    %     zvec(i,1:length(simdata{i}.zvec))=simdata{i}.zvec;
        zvec=simdata{i}.univ.grid.z(:,1,1);
        zvec=zvec-zvec(1); % move such that top of resist is 0

        deprdata=depr_thruz(outparms,elecimg_bin,simdata{i}.acidimg,deprthr);

        zvec2(i,1:length(zvec))=zvec';
        deprmean_tmp(trial_count,1:length(deprdata.mean))=deprdata.mean;
        deprmean2_tmp(trial_count,1:length(deprdata.mean))=deprdata.mean2;
        deprmean3_tmp(trial_count,1:length(deprdata.mean))=deprdata.mean3;

    %     tloop_end=toc(tloop_start);
    %     fprintf('Took %.4f s = %.4f min\n\n',tloop_end,tloop_end/60);
    end
%     zvec2(i,1:length(zvec))=zvec';
    deprmean(i,1:size(deprmean_tmp,2))=mean(deprmean_tmp,1);
    deprmean2(i,1:size(deprmean2_tmp,2))=mean(deprmean2_tmp,1);
    deprmean3(i,1:size(deprmean3_tmp,2))=mean(deprmean3_tmp,1);
    deprstd(i,1:size(deprmean_tmp,2))=std(deprmean_tmp,0,1); % std(X,0,DIM)
    deprstd2(i,1:size(deprmean2_tmp,2))=std(deprmean2_tmp,0,1); % std(X,0,DIM)
    deprstd3(i,1:size(deprmean3_tmp,2))=std(deprmean3_tmp,0,1); % std(X,0,DIM)
end
fprintf('Rcn/Diffn Calculations complete\n\n');
tend=toc(tstart);
fprintf('...%.2f s = %.2f min\n',tend,tend/60);

%% Plotting acid data for publication purposes

close all;

selidx=16;
figure;
imagesc(sum(simdata{selidx}.acidimg,3));
colorbar;
xlim([10 40]);
ylim([1 25]);
caxis([0 10]);
set(gca,'fontsize',30,'linewidth',3.0);
simdata{selidx}.Dose

%% plot the deprotection profiles for publications 
clc;
close all;

thr_tmp=0.5;
outparms=rcndiffndata{16};
tmp=outparms.deprimg(:,:,23:29);
% tmp(outparms.protimg_pre(:,:,23:29)==0)=1;
tmp2=mean(tmp,3);
tmp3=tmp2;
% tmp3=[ones(1,size(tmp3,2));tmp3];
tmp3(tmp3>=thr_tmp)=1;
tmp3(tmp3<thr_tmp)=0;

% figure;
% imagesc(tmp3);colorbar;

% tmp3=imfill(tmp3,'holes');
% tmp3=tmp3(2:end,:);
% tmp3=bwareaopen(tmp3,100);

figure;
imagesc(mean(tmp,3));
colorbar;
% xlim([10 40]);
% ylim([1 25]);
caxis([0 1]);
set(gca,'fontsize',30,'linewidth',3.0);

% figure;
% imagesc(tmp3);colorbar;
% imagesc(edge(tmp3));colorbar;

%% plot thru-z deprotection profiles
% clc;
close all;

ztmp=simdata{1}.univ.grid.z(:,1,1);
 
gcf_pos=[680   223   958   755];
gca_pos=[0.1300    0.1166    0.7750    0.8084];

figure;
% errorbar(ztmp',deprmean(11:end),deprstd(11:end),'-ob')
% errorbar(ztmp-min(ztmp),deprmean2(1,:)',deprstd2(1,:)','-ob','linewidth',3.0,'MarkerSize',6);
plot(ztmp-min(ztmp),deprmean2(10,:)','-ob','linewidth',3.0,'MarkerSize',12);
hold on;
% errorbar(ztmp-min(ztmp),deprmean2(2,:)',deprstd2(2,:)','-ok','linewidth',3.0,'MarkerSize',6);
plot(ztmp-min(ztmp),deprmean2(11,:)','-ok','linewidth',3.0,'MarkerSize',12);
% errorbar(ztmp-min(ztmp),deprmean2(3,:)',deprstd2(2,:)','-or','linewidth',3.0,'MarkerSize',6);
plot(ztmp-min(ztmp),deprmean2(13,:)','-or','linewidth',3.0,'MarkerSize',12);
% errorbar(ztmp-min(ztmp),deprmean2(4,:)',deprstd2(2,:)','-og','linewidth',3.0,'MarkerSize',6);
plot(ztmp-min(ztmp),deprmean2(16,:)','-og','linewidth',3.0,'MarkerSize',12);
xlabel('z (nm)','fontsize',30);
ylabel('Deprotection Level');
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
xlim([0 15]);
ylim([0 0.8]);
box on;


figure;
% errorbar(ztmp',deprmean(11:end),deprstd(11:end),'-ob')
errorbar(ztmp'-min(ztmp),deprmean(2,:)',deprstd(1,:)','-ob','linewidth',3.0,'MarkerSize',6);
hold on;
errorbar(ztmp'-min(ztmp),deprmean(3,:)',deprstd(2,:)','-ok','linewidth',3.0,'MarkerSize',6);
xlabel('z (nm)','fontsize',30);
ylabel('Deprotection Level');
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);
xlim([0 50]);
ylim([0 0.8]);
box on;

%% Create the deprmean matrix from saved rcndiffndata
clc;
close all;
deprthr=0.5;
    deprmeanB=zeros(length(rcndiffndata),50);
    deprstdB=zeros(length(rcndiffndata),50);
    deprmean2B=zeros(length(rcndiffndata),50);
    deprstd2B=zeros(length(rcndiffndata),50);
    deprmean3B=zeros(length(rcndiffndata),50);
    deprstd3B=zeros(length(rcndiffndata),50);
for i = 1:length(rcndiffndata)
    deprmean_tmp=zeros(1,50); % should change if > 50 in the future
    deprmean2_tmp=zeros(1,50);
    deprmean3_tmp=zeros(1,50);
    if ~isempty(rcndiffndata{i})
        deprdata=depr_thruz(rcndiffndata{i},simdata{i}.elecimg_inc,simdata{i}.acidimg,deprthr);
    
    %     zvec2(i,1:length(zvec))=zvec';
    %     deprmean(i,1:length(deprdata.mean))=deprdata.mean;
    %     deprmean2(i,1:length(deprdata.mean))=deprdata.mean2;
    %     deprstd(i,1:length(deprdata.std))=deprdata.mean;

        zvec2(i,1:length(zvec))=zvec';
        deprmean_tmp(1,1:length(deprdata.mean))=deprdata.mean;
        deprmean2_tmp(1,1:length(deprdata.mean))=deprdata.mean2;
        deprmean3_tmp(1,1:length(deprdata.mean))=deprdata.mean3;

    %     tloop_end=toc(tloop_start);
    %     fprintf('Took %.4f s = %.4f min\n\n',tloop_end,tloop_end/60);
    %     zvec2(i,1:length(zvec))=zvec';
        deprmeanB(i,1:size(deprmean_tmp,2))=mean(deprmean_tmp,1);
        deprmean2B(i,1:size(deprmean2_tmp,2))=mean(deprmean2_tmp,1);
        deprmean3B(i,1:size(deprmean2_tmp,2))=mean(deprmean3_tmp,1);
        deprstdB(i,1:size(deprmean_tmp,2))=std(deprmean_tmp,0,1); % std(X,0,DIM)
        deprstd2B(i,1:size(deprmean2_tmp,2))=std(deprmean2_tmp,0,1); % std(X,0,DIM)
        deprstd3B(i,1:size(deprmean2_tmp,2))=std(deprmean3_tmp,0,1); % std(X,0,DIM)
    else
        
    end
end

%% CLEAR Workspace, AND Load existing mat files containing deprmean matrices
clear;
clc;
close all;

filepath='RcnDiffn_Data\';
filename='';
[filename,pathname]=uigetfile([filepath '\*.mat'],'MultiSelect','off');

load([filepath filename]);

%% clear workspace except selected variable names
clc;

clearvars -except fnames pathname simdata

%% clear rtdata;
clear rtdata;

%% Extract thickness loss at the threshold from the deprmean matrix above
clc;
% close all;

%%% 60 eV: dose = 48 is needed.
%%% 80 eV: dose = 24 is needed.
%%% 40 eV: dose = 96 and 100 are needed.

rtlost=[];
rtlost_std=[];
nAcids=[];
rt_offset=[];
depr_thr=0.5;
analysis_entity=deprmean2;

E_choose=[29 40 49 60 70 80 91];
E_choose=[29 40 60];

t_unexp=32;

count=1;
for i = 1:size(analysis_entity,1)
% for i = 86
    ztmp=zvec2(i,:);
    ztmp=[-simparms.zadd_npx:1:ztmp(1)-1 ztmp];
%     ztmp=ztmp-simparms.zadd_npx; % nmpp=1 by default, double check in future
    deprtmp=analysis_entity(i,:);
    deprtmp_std=deprstd2(i,:);
    
    ntrials2=1;
    rt_tmp=[];
    for k = 1:ntrials2
        deprvec=[];
        for j = 1:length(deprtmp)
            deprvec(j)=deprtmp(j)+deprtmp_std(j)*randn;
        end

% % %         algo 1:
% % %         idx=find(deprvec==max(deprvec));
% % %         if ~isempty(idx)
% % %             idx=idx(1);
% % %         end
% % %         if max(deprvec)>depr_thr
% % %             idx3=idx-1+find(deprvec(idx:end)<depr_thr); idx3=idx3(1);
% % %             idx2=idx-1+find(deprvec(idx:end)>depr_thr); idx2=idx2(end);
% % %             p=polyfit(ztmp(min([idx2 idx3]):max([idx2 idx3])),deprvec(min([idx2 idx3]):max([idx2 idx3])),1);
% % %             rt_tmp(k)=(depr_thr-p(2))/p(1);
% % %         else
% % %             rt_tmp(k)=0;
% % %         end

% % %         algo 2:
        idx=find(deprvec(end:-1:1)>depr_thr);
        if ~isempty(idx)
            idx1=length(deprvec)-idx(1)+1;
            p=polyfit([ztmp(idx1) ztmp(idx1+2)],[deprvec(idx1) deprvec(idx1+2)],1);
            rt_tmp(k)=(depr_thr-p(2))/p(1);
        else
            rt_tmp(k)=0;
        end
        
%         Find from beginning of vector as well:
        idx=find(deprvec(1:end)>depr_thr);
        if ~isempty(idx)
            idx1=idx(1);
            idx1=max([idx1-1 1]);
            p=polyfit([ztmp(idx1) ztmp(idx1+1)],[deprvec(idx1) deprvec(idx1+1)],1);
            rt_offset(i,k)=(depr_thr-p(2))/p(1);
            rt_offset(i,k)=max([0 rt_offset(i,k)]);
            rt_tmp(k)=rt_tmp(k)-rt_offset(i,k);
        end
    end
    rtlost(count)=mean(rt_tmp);
    rtlost_std(count)=std(rt_tmp);
    nAcids(count)=sum(simdata{i}.acidimg(:));
    count=count+1;
end

rt=t_unexp-rtlost;
% if ~exist('tmpdata')
%     tmpdata=[];
% end
% tmpdata=[tmpdata rt'];


%%% Thru-Trial data: Break into various doses
% clc;
% close all;

Evec=[];

if ~exist('rtdata')
    rtdata={};
else
%     if isempty(rtdata)
%         rtdata={};
%     else
%         for i = 1:length(rtdata)
%             fieldnames=fields(rtdata{i});
%             for j = 1:length(fieldnames)
%                 rtdata{i}=setfield(rtdata{i},fieldnames{j},[]);
%             end
%         end
%     end
end

dosevec=[];
for i = 1:length(simdata)
    Evec(i)=simdata{i}.E;
    dosevec(i)=simdata{i}.Dose/prod(simdata{i}.univ.px_nm(2:3));
end

dose_scale=[0.45 0.53 0.83 1]; % [29 eV, 49 eV, 70 eV, 91 eV]
% dose_scale=[0.45 0.53 1]; % [29 eV, 49 eV, 70 eV, 91 eV]

for i = 1:length(E_choose)
    idx=find(Evec==E_choose(i));
    dose_tmp=dosevec(idx);
    nAcids_tmp=nAcids(idx);
    rt_tmp=t_unexp-rtlost(idx);
    std_tmp=rtlost_std(idx);
    
    if ~isempty(dose_tmp)
        dose_unique=unique(dose_tmp);
    else
        error('ERROR: dose_tmp is empty');
    end
    if depr_thr==0.3
        rowsel=1;
    else
        rowsel=2;
    end
    rt_mean=[];
    rt_std=[];
    for j = 1:length(dose_unique)
        idx2=find(dose_tmp==dose_unique(j));
        rt_mean(j)=mean(rt_tmp(idx2));
%         rt_std(j)=std(rt_tmp(idx2));
        rt_std(j)=sqrt((std(rt_tmp(idx2)))^2+(mean(std_tmp(idx2)))^2);
    end
    
    rtdata{i}.nAcids=nAcids_tmp;
    rtdata{i}.dose=dose_tmp;
    rtdata{i}.rt=rt_tmp;
    rtdata{i}.dosemean=dose_unique;
%     rtdata{i}.dosemean=dose_unique*dose_scale(i);
    rtdata{i}.E=E_choose(i);
    
    rtdata{i}.rtmean(rowsel,1:length(rt_mean))=rt_mean;
    rtdata{i}.rtstd(rowsel,1:length(rt_std))=rt_std;
end

figure;
plot(Evec);title('E');

% rtdata{2}
% figure;
% plot(dosevec);title('dose');

%% Plot the thru-trial data stored in "rtdata" above:
% clc;
close all;

%%% first divide up the rt data into rt_0p3 and rt_0p5 for the 2 thresholds
for i = 1:length(rtdata)
    rtdata{i}.rtmean_0p3=rtdata{i}.rtmean(1,:);
    rtdata{i}.rtstd_0p3=rtdata{i}.rtstd(1,:);
    
    rtdata{i}.rtmean_0p5=rtdata{i}.rtmean(2,:);
    rtdata{i}.rtstd_0p5=rtdata{i}.rtstd(2,:);
end
p_0p3=[];
p_0p5=[];

gca_pos=[0.1300    0.1100    0.7750    0.8150];
gcf_pos=[680   366   799   612];
pltstyles={'-db','-ok','-or','-og','--ob','--ok','--or','--og'};
pltstyles2={'--b','--k','--r','--g','--ob','--ok','--or','--og'};

% figure;hold on;xlabel('Dose (e^-/nm^2)');ylabel('Thickness (nm)');
% set(gca,'fontsize',30,'linewidth',3.0,'XScale','log');
% set(gcf,'Position',gcf_pos);
% set(gca,'Position',gca_pos);
% 
% for i = 1:length(rtdata)
%     dosevec=rtdata{i}.dosemean;
%     rtvec=rtdata{i}.rtmean_0p3;
%     rt_lb=38;    rt_ub=44;
%     p=polyfit(log10(dosevec(rtvec<rt_ub&rtvec>=rt_lb)),rtvec(rtvec<rt_ub&rtvec>=rt_lb),1);
%     p_0p3(i)=p(1);
%     
%     errorbar(rtdata{i}.dosemean,rtdata{i}.rtmean_0p3,rtdata{i}.rtstd_0p3,pltstyles{i},'linewidth',3.0,'markersize',12);
%     drawnow;
% end

%%%% threshold=0.5
fig1=figure;hold on;xlabel('Dose (e^-/nm^2)');ylabel('Thickness (nm)');
set(gca,'fontsize',30,'linewidth',3.0,'XScale','log');
box on;
set(gcf,'Position',gcf_pos);
set(gca,'Position',gca_pos);

% fig2=figure;hold on;xlabel('Dose (e^-/nm^2)');ylabel('# of Acids');
% set(gca,'fontsize',30,'linewidth',3.0,'XScale','linear');
% box on;
% set(gcf,'Position',gcf_pos);
% set(gca,'Position',gca_pos);

for i = 1:length(rtdata)
% for i = 1
    dosevec=rtdata{i}.dosemean;
    rtvec=rtdata{i}.rtmean_0p5;
    idx=find(rtvec==t_unexp);
    if ~isempty(idx)
        idx=idx(end);
    else
        idx=1;
    end
    dosevec=dosevec(idx:end);
    rtvec=rtvec(idx:end);
    
    if E_choose(i)==40
        rt_lb=40-(50-t_unexp);    rt_ub=49.95-(50-t_unexp); % 40 eV
        rt_lb=41;    rt_ub=50; % 40 eV; thr=0.5
        rt_lb=41;    rt_ub=49.96; % 40 eV; thr=0.5; 4 Rcn/DiffnTrials; exclude last few points
        rt_lb=40;    rt_ub=50; % 40 eV; thr=0.5; 4 Rcn/DiffnTrials simrun; include all points
%         rt_lb=42;    rt_ub=50; % 40 eV; thr=0.55;
%         rtvec=rtvec(dosevec<=256);      % exclude last few points
%         dosevec=dosevec(dosevec<=256);  % exclude last few points
    end
    if E_choose(i)==60
        rt_lb=41.6-(50-t_unexp);    rt_ub=49.96-(50-t_unexp); % 60 eV
        rt_lb=41.6-(50-t_unexp);    rt_ub=50-(50-t_unexp); % 60 eV
        rt_lb=38-(50-t_unexp);    rt_ub=50-(50-t_unexp); % 60 eV
        rt_lb=40;    rt_ub=50; % 60 eV; thr=0.5
        rt_lb=40;    rt_ub=50; % 60 eV; thr=0.5; 4 Rcn/DiffnTrials; exclude last few points
        rt_lb=38;    rt_ub=50; % 40 eV; thr=0.5; 4 Rcn/DiffnTrials simrun; include all points
%         rt_lb=42;    rt_ub=50; % 60 eV; thr=0.55
%         rtvec=rtvec(dosevec<=128);      % exclude last few points
%         dosevec=dosevec(dosevec<=128);  % exclude last few points
    end
    if E_choose(i)==80
        rt_lb=40-(50-t_unexp);    rt_ub=49.96-(50-t_unexp); % 80 eV
        rt_lb=42.5;    rt_ub=50; % 80 eV; thr = 0.5
        rt_lb=42;    rt_ub=49.2; % 80 eV; thr = 0.5; 4 Rcn/DiffnTrials simrun
        rt_lb=38;    rt_ub=49.7; % 80 eV; thr=0.5; 4 Rcn/DiffnTrials simrun; exclude last few points
        rt_lb=41;    rt_ub=50; % 40 eV; thr=0.5; 4 Rcn/DiffnTrials simrun; include all points
%         rtvec=rtvec(dosevec<=128);      % exclude last few points
%         dosevec=dosevec(dosevec<=128);  % exclude last few points
    end
    if E_choose(i)==29
        rt_lb=0; rt_ub=31;
        rt_lb=0; rt_ub=32;
%         rtvec=rtvec(dosevec<=8);      % exclude last few points
%         dosevec=dosevec(dosevec<=8);  % exclude last few points
    end
    if E_choose(i)==49
        rt_lb=0; rt_ub=31;
        rt_lb=0; rt_ub=32;
    end
    if E_choose(i)==91
        rt_lb=0; rt_ub=31;
        rt_lb=0; rt_ub=32;
%         rtvec=rtvec(dosevec<=6);      % exclude last few points
%         dosevec=dosevec(dosevec<=6);  % exclude last few points
    end
    dosefit=dosevec(rtvec<=rt_ub&rtvec>=rt_lb);
    rtfit=rtvec(rtvec<=rt_ub&rtvec>=rt_lb);
%     rtfit=rtfit(dosefit==32|dosefit==64|dosefit==128|dosefit==256);
%     dosefit=dosefit(dosefit==32|dosefit==64|dosefit==128|dosefit==256);
    p=polyfit(log10(dosefit),rtfit,1);
    p_0p5(i)=p(1);
    
    figure(fig1);
%     errorbar(rtdata{i}.dosemean,rtdata{i}.rtmean_0p5,rtdata{i}.rtstd_0p5,pltstyles{i},'linewidth',3.0,'markersize',12);
    errorbar(rtdata{i}.dosemean.*(rtdata{i}.E-0),rtdata{i}.rtmean_0p5,rtdata{i}.rtstd_0p5,pltstyles{i},'linewidth',3.0,'markersize',12);
%     plot(dosefit,polyval(p,log10(dosefit)),pltstyles2{i},'linewidth',5.0);
%     xlim([0.1 100]);
%     ylim([17 32]);
    drawnow;
    
%     figure(fig2);
%     plot(rtdata{i}.dose,rtdata{i}.nAcids,pltstyles{i},'linewidth',3.0,'markersize',12);
% %     plot(rtdata{i}.dose./rtdata{i}.dose(1),rtdata{i}.nAcids./rtdata{i}.nAcids(1),pltstyles{i},'linewidth',3.0,'markersize',12);
%     drawnow;
end
% figure(fig2); plot([1:16],1:16,'--k','linewidth',3.0);

[p_0p3;p_0p5]

%% Generalized data structure to hold dose, thickness, slope values
% clc;
close all;

contrcurve_data=gensimdata(dose,t_unexp-rtlost,rtlost_std);
% contrcurve_data.slope=polyfit(log10(contrcurve_data.dose(4:end)./contrcurve_data.dose(4)),contrcurve_data.rt(4:end),1);
contrcurve_data.slope=polyfit(log10(contrcurve_data.dose(1:end)),contrcurve_data.rt(1:end),1);

figure;plot(contrcurve_data.dose,contrcurve_data.rt,'-o');
xlim([0.9*min(contrcurve_data.dose) 1.1*max(contrcurve_data.dose)]);
set(gca,'XScale','log');
contrcurve_data.slope

%% create the data structures that hold dose, thickness values
sim29={};
sim49={};
sim91={};
sim29.dose=dose(1,1:5);
sim29.rt=t_unexp-rtlost(1,1:5);
sim29.rt_eb=rtlost_std(1,1:5);
sim29.slope=polyfit(log10(sim29.dose),sim29.rt,1)
% try
%     sim49.dose=dose(1,6:10);
%     sim49.rt=t_unexp-rtlost(1,6:10);
%     sim49.rt_eb=rtlost_std(1,6:10);
%     sim49.slope=polyfit(log10(sim49.dose),sim49.rt,1);
% catch
% %     fprintf('error in sim49');
% end

try
    sim91.dose=dose(1,6:9);
    sim91.rt=t_unexp-rtlost(1,6:9);
    sim91.rt_eb=rtlost_std(1,6:9);
    sim91.slope=polyfit(log10(sim91.dose),sim91.rt,1)
catch
end

%%% find where rtlost==tunexposed
% simulation:
rt_sim=t_unexp-rtlost(1:end);
dose_sim=dose(1:end);

% experiment: 29 eV
try
    data29.dose2=unique(data29.dose);
    for i = 1:length(data29.dose2)
        data29.rt2_mean(i,1)=mean(data29.rt(data29.dose==data29.dose2(i)));
        data29.rt2_std(i,1)=std(data29.rt(data29.dose==data29.dose2(i)));
    end
catch end

try
    data49.dose2=unique(data49.dose);
    for i = 1:length(data49.dose2)
        data49.rt2_mean(i,1)=mean(data49.rt(data49.dose==data49.dose2(i)));
        data49.rt2_std(i,1)=std(data49.rt(data49.dose==data49.dose2(i)));
    end
catch end

try
    data91.dose2=unique(data91.dose);
    for i = 1:length(data91.dose2)
        data91.rt2_mean(i,1)=mean(data91.rt(data91.dose==data91.dose2(i)));
        data91.rt2_std(i,1)=std(data91.rt(data91.dose==data91.dose2(i)));
    end
catch end

%% Analyze only the simulation results
% clc;
% close all;

% sim29.rt=sim29.rt(sim29.dose~=25);
% sim29.dose=sim29.dose(sim29.dose~=25);
% sim29.rt_eb=sim29.rt_eb(sim29.dose~=25);
% 
% sim49.rt=sim49.rt(sim49.dose~=25);
% sim49.dose=sim49.dose(sim49.dose~=25);
% sim49.rt_eb=sim49.rt_eb(sim49.dose~=25);
% 
% sim91.rt=sim91.rt(sim91.dose~=25);
% sim91.dose=sim91.dose(sim91.dose~=25);
% sim91.rt_eb=sim91.rt_eb(sim91.dose~=25);
% 
% sim29_2.rt=sim29_2.rt(sim29_2.dose~=25);
% sim29_2.dose=sim29_2.dose(sim29_2.dose~=25);
% sim29_2.rt_eb=sim29_2.rt_eb(sim29_2.dose~=25);
% 
% sim49_2.rt=sim49_2.rt(sim49_2.dose~=25);
% sim49_2.dose=sim49_2.dose(sim49_2.dose~=25);
% sim49_2.rt_eb=sim49_2.rt_eb(sim49_2.dose~=25);
% 
% sim91_2.rt=sim91_2.rt(sim91_2.dose~=25);
% sim91_2.dose=sim91_2.dose(sim91_2.dose~=25);
% sim91_2.rt_eb=sim91_2.rt_eb(sim91_2.dose~=25);

gcf_pos=[680   236   926   742];
gca_pos=[  0.1300    0.1612    0.7750    0.7638];

p1=polyfit(log10(sim29.dose),sim29.rt,1);
try
    p2=polyfit(log10(sim91.dose),sim91.rt,1);
catch
end

figure;
set(gcf,'Position',gcf_pos);
set(gca,'Position',gca_pos);
errorbar(1.*sim29.dose,sim29.rt,sim29.rt_eb,'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
% plot(sim29.dose,polyval(p1,log10(sim29.dose)),'k','linewidth',3.0);
% errorbar(1.*sim49.dose,sim49.rt,sim49.rt_eb,'-ok','linewidth',3.0,'MarkerSize',12);
try
    errorbar(sim91.dose,sim91.rt,sim91.rt_eb,'-or','linewidth',3.0,'MarkerSize',12);
catch
end
% plot(sim91.dose,polyval(p2,log10(sim91.dose)),'k','linewidth',3.0);
xlabel('Dose (e^-/nm^2)','fontsize',30);
ylabel('Thicknss (nm)','fontsize',30);
set(gca,'fontsize',30,'linewidth',3.0,'XScale','log');
xlim([0.5 20]);
ylim([38 50]);
legend('29 eV','49 eV','91 eV');
legend('29 eV','91 eV');
title(sprintf('thr = %.2f; SLopes = [%.2f %.2f]',depr_thr,sim29.slope(1),sim91.slope(1)));

% figure;
% set(gcf,'Position',gcf_pos);
% set(gca,'Position',gca_pos);
% errorbar(0.33.*sim29.dose,sim29.rt,sim29.rt_eb,'-ob','linewidth',3.0,'MarkerSize',12);
% hold on;
% % errorbar(49/91.*sim49.dose,sim49.rt,sim49.rt_eb,'-ok','linewidth',3.0,'MarkerSize',12);
% errorbar(sim91.dose,sim91.rt,sim91.rt_eb,'-or','linewidth',3.0,'MarkerSize',12);
% xlabel('Dose (e^-/nm^2)','fontsize',30);
% ylabel('Thicknss (nm)','fontsize',30);
% set(gca,'fontsize',30,'linewidth',3.0,'XScale','linear');
% % xlim([0.7 300]);
% % ylim([40 52]);
% legend('29 eV','49 eV','91 eV');

%% Total # of acids in the LEEM simulations
% clc;
% close all;

nAcids=[];
zval=[];
dose=[];
E=[];
count=1;
% for i = 1:length(simdata)
% for i = 1:8 % 29 eV
for i = 9:16 % 40 eV
% for i = 15:20 % 60 eV
% for i = 21:26 % 80 eV
% for i = 27:31 % 91 eV
    if simdata{i}.zval~=0
%     nAcids(i)=sum(simdata{i}.acidimg(:));
        nAcids(count)=(simdata{i}.nAcids);
        zval(count)=simdata{i}.zval;
        dose(count)=simdata{i}.Dose;
        E(count)=simdata{i}.E;
        count=count+1;
    end
end
[dose,idx]=sort(dose);
nAcids=nAcids(idx);

gcf_pos=[680   359   799   619];
gca_pos=[0.1300    0.1100    0.7750    0.8150];
% figure;
figure(fig1);
plot(dose,nAcids,'--ok','linewidth',3.0,'markersize',12); % total number
% plot(dose,nAcids./(36.*dose),'--ok','linewidth',3.0,'markersize',12); % per incident electron
% ylim([0 300]);xlim([0 32]);
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'XScale','linear');
set(gcf,'Position',gcf_pos);
set(gca,'Position',gca_pos);
box on;

figure;
plot(nAcids,'-ob','markersize',12);
hold on;
plot([1 length(nAcids)],[1 1].*mean(nAcids),'--b','linewidth',3.0);
plot([1 length(nAcids)],[1 1].*(mean(nAcids)+std(nAcids)),'--k','linewidth',3.0);
plot([1 length(nAcids)],[1 1].*(mean(nAcids)-std(nAcids)),'--k','linewidth',3.0);
plot([1 length(nAcids)],[1 1].*(mean(nAcids)+2*std(nAcids)),'--r','linewidth',3.0);
plot([1 length(nAcids)],[1 1].*(mean(nAcids)-2*std(nAcids)),'--r','linewidth',3.0);
plot([1 length(nAcids)],[1 1].*(mean(nAcids)+3*std(nAcids)),'--g','linewidth',3.0);
plot([1 length(nAcids)],[1 1].*(mean(nAcids)-3*std(nAcids)),'--g','linewidth',3.0);
title(sprintf('Mean = %.2f; Std = %.2f',mean(nAcids),std(nAcids)));
figure;
plot(zval);

%% Saturation tests
clc;
close all;

Evec_acidsim=[];
dosevec_acidsim=[];
nacids_sat=[];
nacids_unsat=[];

nacids_total=[];
nacids_znorm=[];
nacids_perelec=[];
depth=[];
acidsimdata={};

E_choose=[29 49 91];

%%%% get data, put in global vectors 
for i = 1:length(simdata)
    
    Evec_acidsim(i)=simdata{i}.E;
    dosevec_acidsim(i)=simdata{i}.Dose;
    nacids_sat(i)=simdata{i}.nacids_thrutrial;
%     nacids_unsat(i)=simdata{i}.nacids_unsat_thrutrial;
    
end

%%%% isolate data for each energy specified in Evec_acidsim
Evec_unique=unique(Evec_acidsim);

for i = 1:length(E_choose)
    idx=find(Evec_acidsim==E_choose(i));
    acidsimdata{i}.E=E_choose(i);
    acidsimdata{i}.nacids_sat=nacids_sat(idx);
    acidsimdata{i}.nacids_unsat=nacids_unsat(idx);
    acidsimdata{i}.nacids_ratio=nacids_unsat(idx)./nacids_sat(idx);
    acidsimdata{i}.dose=dosevec_acidsim(idx);
end

%% Auxiliary: Write contents of fnames into a file with numbering
clc;
close all;

fid=fopen('Fnames_List.dat','w');

for i = 1:length(fnames)
    fprintf(fid,'%d \t %s\n',i,fnames{i});
end

fclose(fid);

%% # of Acids thru-z in the LEEM simulations
% clc;
close all;

Acidstats=load('QE_Stats_low10PAG\AcidStats_thruE.mat');
Acidstats2=load('QE_low10PAG_ElecOnTop\AcidStats_thruE.mat');

nAcids_thruz=[];
Evec_acidsim=[];
dosevec_acidsim=[];
nacids_total=[];
nacids_znorm=[];
nacids_perelec=[];
totalacids_sat=[];
nacids_perE=[];
depth=[];
acidsimdata={};

E_choose=[29 49 70 91];

dose_scale=[0.45 0.53 0.83 1]; % [29 eV, 49 eV, 70 eV, 91 eV]

%%%% get data, put in global vectors 
for i = 1:length(simdata)
    acidimg=simdata{i}.acidimg;
    acidimg_zsum=reshape(sum(acidimg,1),50,50);
    idx=find(acidimg_zsum>0);
    [xidx,yidx]=ind2sub(size(acidimg_zsum),idx);
    slice_area=(max(xidx)-min(xidx))*(max(yidx)-min(yidx));
        
%     acid_pernm2=[];
%     for j = 1:size(acidimg,1)
%         acidimg_slice=reshape(acidimg(j,:,:),50,50);
%         binimg=double((acidimg_slice>1));
%         idx=find(acidimg_slice>0);
%         [xidx,yidx]=ind2sub(size(acidimg_slice),idx);
% %         slice_area=(max(xidx)-min(xidx))*(max(yidx)-min(yidx));
%         acid_pernm2(j)=sum(acidimg_slice)/slice_area;
%     end
    
    Evec_acidsim(i)=simdata{i}.E;
    dosevec_acidsim(i)=simdata{i}.Dose;
    nacids_total(i)=sum(acidimg(:));
    nacids_perelec(i)=nacids_total(i)/(dosevec_acidsim(i)*36);
    nacids_perE(i)=nacids_total(i)/(simdata{i}.E);
    nacids_pere_perE(i)=nacids_total(i)/(simdata{i}.E.*dosevec_acidsim(i)*36);
    
    totalacids_sat(i)=simdata{i}.nacids_unsat_total;
    for j = 1:size(acidimg,1)
        nAcids_thruz(i,j)=sum(sum(sum(acidimg(j,:,:))));
    end
    idx=find(nAcids_thruz(i,:)==0);
    depth(i)=idx(1);
    nacids_znorm(i)=nacids_perelec(i)/depth(i);
end

%%%% isolate data for each energy specified in Evec_acidsim
Evec_unique=unique(Evec_acidsim);

for i = 1:length(E_choose)
    idx=find(Evec_acidsim==E_choose(i));
    acidsimdata{i}.nAcids=nAcids_thruz(idx,:);
    acidsimdata{i}.dose=dosevec_acidsim(idx);
    acidsimdata{i}.eVnm2=acidsimdata{i}.dose.*E_choose(i);
    
    acidsimdata{i}.totalacids=nacids_total(idx);
    acidsimdata{i}.totalacids_sat=totalacids_sat(idx);
    acidsimdata{i}.acids_znorm=nacids_znorm(idx);
    acidsimdata{i}.nacids_perelec=nacids_perelec(idx);
    acidsimdata{i}.nacids_perE=nacids_perE(idx);
    acidsimdata{i}.nacids_pere_perE=nacids_pere_perE(idx);
    
    acidsimdata{i}.depth=depth(idx);
    acidsimdata{i}.dose_scale=dose_scale(i);
    
    acidsimdata{i}.dose_unique=unique(acidsimdata{i}.dose);
    acidsimdata{i}.eVnm2_unique=unique(acidsimdata{i}.eVnm2);
    
    for j = 1:length(acidsimdata{i}.dose_unique)
        idx2=find(acidsimdata{i}.dose==acidsimdata{i}.dose_unique(j));
        
        acidsimdata{i}.nacids_mean(j)=mean(acidsimdata{i}.totalacids(idx2));
        acidsimdata{i}.nacids_std(j)=std(acidsimdata{i}.totalacids(idx2));
        
        acidsimdata{i}.acdistat_mean(j)=mean(acidsimdata{i}.totalacids_sat(idx2));
        acidsimdata{i}.acdistat_std(j)=std(acidsimdata{i}.totalacids_sat(idx2));
        
        acidsimdata{i}.yield_mean(j)=mean(acidsimdata{i}.totalacids(idx2)./acidsimdata{i}.totalacids_sat(idx2));
        acidsimdata{i}.yield_std(j)=std(acidsimdata{i}.totalacids(idx2)./acidsimdata{i}.totalacids_sat(idx2));
        
        acidsimdata{i}.yield_perE_mean(j)=mean(acidsimdata{i}.totalacids(idx2)./E_choose(i));
        acidsimdata{i}.yield_perE_std(j)=std(acidsimdata{i}.totalacids(idx2)./E_choose(i));

        acidsimdata{i}.nacids_pe_mean(j)=mean(acidsimdata{i}.nacids_perelec(idx2));
        acidsimdata{i}.nacids_pe_std(j)=std(acidsimdata{i}.nacids_perelec(idx2));
    end
end

dosevec_global=[];
for i = 1:length(acidsimdata)
    dosevec_global=[dosevec_global unique(acidsimdata{i}.dose)];
end

dosevec_global=unique(dosevec_global);
acids_perelec_thruE=[];
acids_perE_mean=[];
acids_perE_std=[];
acids_znorm_thruE=[];
totalacids_thruE=[];
acids_thruE={};
for i = 1:length(dosevec_global)
    for j = 1:length(acidsimdata)
        idx=find(acidsimdata{j}.dose==dosevec_global(i));
        if ~isempty(idx)
            acids_perelec_thruE(i,j)=mean(acidsimdata{j}.nacids_perelec(idx));
            acids_perE_mean(j,i)=mean(acidsimdata{j}.nacids_perE(idx));
            acids_perE_std(j,i)=std(acidsimdata{j}.nacids_perE(idx));
            
            acids_pere_perE_mean(j,i)=mean(acidsimdata{j}.nacids_pere_perE(idx));
            acids_pere_perE_std(j,i)=std(acidsimdata{j}.nacids_pere_perE(idx));
            
            totalacids_thruE(i,j)=mean(acidsimdata{j}.totalacids(idx));
            acids_znorm_thruE(i,j)=mean(acidsimdata{j}.acids_znorm(idx));
            acids_thruE{i}.zdist(j,:)=mean(acidsimdata{j}.nAcids(idx,:),1);
            acids_thruE{i}.dose(j)=mean(acidsimdata{j}.dose(idx));
        else
            acids_perelec_thruE(i,j)=NaN;
            totalacids_thruE(i,j)=NaN;
            acids_znorm_thruE(i,j)=NaN;
            acids_thruE{i}.zdist(j,:)=NaN.*ones(1,50);
            acids_thruE{i}.dose(j)=NaN;
        end
    end
end

%%%% plot the acids/electron vs. dose for all electron energies
gcf_pos=[322   300   918   678];
gca_pos=[0.1300    0.1195    0.7750    0.8055];

pltstyle={'-db','-dk','-dr','-dg','-ob','-ok','-or'};

%%%% Nacids vs. e/nm2
figure;hold on;box on;
xlabel('Dose (e^-/nm^2)');ylabel('Total Acids');
set(gca,'XScale','log','fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);set(gca,'Position',gca_pos);
for i = 1:length(acidsimdata)
%     plot(acidsimdata{i}.dose,acidsimdata{i}.totalacids,pltstyle{i},'linewidth',3.0,'markersize',12);
    errorbar(acidsimdata{i}.dose_unique,acidsimdata{i}.nacids_mean,acidsimdata{i}.nacids_std,pltstyle{i},'linewidth',3.0,'markersize',12);
    drawnow;
end
xlim([0.1 100]);

%%%% Nacids vs. eV/nm2
figure;hold on;box on;
xlabel('Dose (eV/nm^2)');ylabel('Total Acids');
set(gca,'XScale','log','fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);set(gca,'Position',gca_pos);
for i = 1:length(acidsimdata)
%     plot(acidsimdata{i}.dose,acidsimdata{i}.totalacids,pltstyle{i},'linewidth',3.0,'markersize',12);
    errorbar(acidsimdata{i}.eVnm2_unique,acidsimdata{i}.nacids_mean,acidsimdata{i}.nacids_std,pltstyle{i},'linewidth',3.0,'markersize',12);
    drawnow;
end
% xlim([0.1 100]);

%%%% Nacids per eV vs. e/nm2
figure;hold on;box on;
xlabel('Dose (e^-/nm^2)');ylabel('Acids per eV');
set(gca,'XScale','log','fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);set(gca,'Position',gca_pos);
for i = 1:length(acidsimdata)
%     plot(acidsimdata{i}.dose,acidsimdata{i}.totalacids,pltstyle{i},'linewidth',3.0,'markersize',12);
    errorbar(acidsimdata{i}.dose_unique,acidsimdata{i}.yield_perE_mean,acidsimdata{i}.yield_perE_std,pltstyle{i},'linewidth',3.0,'markersize',12);
    drawnow;
end
% xlim([0.1 100]);

%%%% Nacids per electron per eV vs. e/nm2
figure;hold on;box on;
xlabel('Energy (eV)');ylabel('Acids per e^- eV');
set(gca,'XScale','linear','fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);set(gca,'Position',gca_pos);
for i = 1:size(acids_pere_perE_mean,2)
%     plot(acidsimdata{i}.dose,acidsimdata{i}.totalacids,pltstyle{i},'linewidth',3.0,'markersize',12);
    errorbar(E_choose',acids_pere_perE_mean(:,i),acids_pere_perE_std(:,i),pltstyle{i},'linewidth',3.0,'markersize',12);
    drawnow;
end
% xlim([0.1 100]);

% % figure;hold on;box on;
% % xlabel('Scaled Dose (e^-/nm^2)');ylabel('Total Acids');
% % set(gca,'XScale','log','fontsize',30,'linewidth',3.0);
% % set(gcf,'Position',gcf_pos);set(gca,'Position',gca_pos);
% % for i = 1:length(acidsimdata)
% % %     plot(acidsimdata{i}.dose.*acidsimdata{i}.dose_scale,acidsimdata{i}.totalacids,pltstyle{i},'linewidth',3.0,'markersize',12);
% %     errorbar(acidsimdata{i}.dose_unique.*acidsimdata{i}.dose_scale,acidsimdata{i}.nacids_mean,acidsimdata{i}.nacids_std,pltstyle{i},'linewidth',3.0,'markersize',12);
% %     drawnow;
% % end
% % xlim([0.1 100]);

figure;hold on;box on;
xlabel('Dose (e^-/nm^2)');ylabel('Acids/Electron');
set(gca,'XScale','log','fontsize',30,'linewidth',3.0);
set(gcf,'Position',gcf_pos);set(gca,'Position',gca_pos);
for i = 1:length(acidsimdata)
%     plot(acidsimdata{i}.dose,acidsimdata{i}.totalacids,pltstyle{i},'linewidth',3.0,'markersize',12);
    errorbar(acidsimdata{i}.dose_unique,acidsimdata{i}.nacids_pe_mean,acidsimdata{i}.nacids_pe_std,pltstyle{i},'linewidth',3.0,'markersize',12);
    drawnow;
end
xlim([0.1 100]);

% % figure;hold on;box on;
% % xlabel('Dose (e^-/nm^2)');ylabel('Acid Yield Ratio');
% % set(gca,'XScale','log','fontsize',30,'linewidth',3.0);
% % set(gcf,'Position',gcf_pos);set(gca,'Position',gca_pos);
% % for i = 1:length(acidsimdata)
% % %     plot(acidsimdata{i}.dose.*acidsimdata{i}.dose_scale,acidsimdata{i}.totalacids./acidsimdata{i}.totalacids_sat,pltstyle{i},'linewidth',3.0,'markersize',12);
% %     errorbar(acidsimdata{i}.dose_unique.*acidsimdata{i}.dose_scale,acidsimdata{i}.yield_mean,acidsimdata{i}.yield_std,pltstyle{i},'linewidth',3.0,'markersize',12);
% % %     plot(acidsimdata{i}.dose_unique.*acidsimdata{i}.dose_scale,acidsimdata{i}.nacids_mean./acidsimdata{i}.acdistat_mean,pltstyle{i},'linewidth',3.0,'markersize',12);
% %     drawnow;
% % end
% % xlim([0.1 100]);
%% rest of the plots from previous cell
figure;
plot(E_choose',acids_perelec_thruE','-o');
hold on;
plot(Acidstats.E',Acidstats.meanAcids','--b','linewidth',3.0);
xlabel('E (eV)');ylabel('Acids/electron');

figure;
plot(E_choose',totalacids_thruE','-o');
hold on;
% plot(Acidstats.E',Acidstats.meanAcids','--b','linewidth',3.0);
xlabel('E (eV)');ylabel('Total Acids');

figure;
plot(E_choose',acids_perelec_thruE'./repmat(acids_perelec_thruE(:,1),1,size(acids_perelec_thruE,2))','-o');
hold on;
plot(E_choose',acids_perelec_thruE(1,:)./acids_perelec_thruE(1,1),'-b','linewidth',3.0);
plot(E_choose',acids_perelec_thruE(end,:)./acids_perelec_thruE(end,1),'-k','linewidth',3.0);
interpval=interp1(Acidstats.E,Acidstats.meanAcids,29);
plot(Acidstats.E',Acidstats.meanAcids'./interpval,'--b','linewidth',3.0);
xlabel('E (eV)');ylabel('Acids/electron');

figure;
plot(E_choose',acids_znorm_thruE'./repmat(acids_znorm_thruE(:,1),1,size(acids_znorm_thruE,2))','-o');
hold on;
plot(E_choose',acids_znorm_thruE(1,:)./acids_znorm_thruE(1,1),'-b','linewidth',3.0);
plot(E_choose',acids_znorm_thruE(end,:)./acids_znorm_thruE(end,1),'-k','linewidth',3.0);
title('Z-Norm Nacids/depth');

figure;
plot(Acidstats.E',Acidstats.meanAcids','-ob','linewidth',3.0);
hold on;
plot(Acidstats2.E',Acidstats2.meanAcids','-ok','linewidth',3.0);

% figure;
% plot(acidsimdata{1}.nAcids(3,:),'-b','linewidth',3.0);
% hold on;
% plot(acidsimdata{2}.nAcids(4,:),'-k','linewidth',3.0);
% plot(acidsimdata{3}.nAcids(4,:),'-r','linewidth',3.0);

% % pltstyles={'-b','-k','-r'};
% % for i = 1:length(acids_thruE)
% %     figure;hold on;
% %     for j = 1:size(acids_thruE{i}.zdist,1)
% %         plot(acids_thruE{i}.zdist(j,:),pltstyles{j},'linewidth',3.0);
% %     end
% %     title(sprintf('Dose = %.2f',acids_thruE{i}.dose(end)));
% %     ylim([0 70]);
% %     xlim([0 20]);
% %     drawnow;
% % end

%% load files into tmp1,tmp2,tmp3 etc... then look at act_global
% clc;
close all;

act_acid=0;
tmp1=load([pathname fnames{3}]);
for i = 1:length(tmp1.act_global)
    if strcmp(tmp1.act_global{i},'acid')==1
        act_acid=act_acid+1;
    end
end
act_acid
