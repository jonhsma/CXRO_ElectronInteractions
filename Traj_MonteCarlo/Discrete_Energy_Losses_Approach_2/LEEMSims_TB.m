%% LEEM simulations analysis
clear;
clc;
close all;

simdata={};

[fname_tmp,pathname]=uigetfile(['LEEM_SimResults_2\*.mat'],'MultiSelect','on');

if ~iscell(fname_tmp)
    fnames{1}=(fname_tmp);
else
    fnames=fname_tmp;
end

elecimg_bin=[];
for i = 1:length(fnames)
    tmp=load([pathname fnames{i}]);
    simdata{i}.acidimg=tmp.acidimg;
    simdata{i}.univ=tmp.univ;
    simdata{i}.elecimg_inc=tmp.elecimg_inc;
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

% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=15eV_Dose=1_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=15eV_Dose=2_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=15eV_Dose=3_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=15eV_Dose=5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=15eV_Dose=10_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=15eV_Dose=20_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');

% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=20eV_Dose=0.25_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=20eV_Dose=0.5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=20eV_Dose=1_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=20eV_Dose=2_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=20eV_Dose=3_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=20eV_Dose=5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=20eV_Dose=10_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=20eV_Dose=20_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');

% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=0.5_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=1_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=2_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=3_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=5_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=10_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=20_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=50_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');

% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=0.25_RCNRAD=2nm_rhoPAG=0.2_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=0.5_RCNRAD=2nm_rhoPAG=0.2_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=1_RCNRAD=2nm_rhoPAG=0.2_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=2_RCNRAD=2nm_rhoPAG=0.2_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=3_RCNRAD=2nm_rhoPAG=0.2_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=5_RCNRAD=2nm_rhoPAG=0.2_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=10_RCNRAD=2nm_rhoPAG=0.2_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');

% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=0.25_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=0.5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=1_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=2_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=3_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=10_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');

% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=20_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=10_RCNRAD=2nm_4nmX4nmEBeam_Univ=30nmX30nmX30nm_px=2nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=50_RCNRAD=2nm_4nmX4nmEBeam_Univ=30nmX30nmX30nm_px=2nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=50_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=40eV_Dose=100_RCNRAD=2nm_4nmX4nmEBeam_Univ=30nmX30nmX30nm_px=2nm.mat'); % THis is weird with thickness loss, due to pixel size discrepancy

% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.1_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm_randPAGs.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.25_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm_randPAGs.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm_randPAGs.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=1_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm_randPAGs.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=2_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm_randPAGs.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=3_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm_randPAGs.mat');

% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.1_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.25_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=1_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=2_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=3_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=5_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=10_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=20_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=50_RCNRAD=2nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');

% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.25_RCNRAD=4nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.5_RCNRAD=4nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=1_RCNRAD=4nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=2_RCNRAD=4nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=3_RCNRAD=4nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=5_RCNRAD=4nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=10_RCNRAD=4nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=20_RCNRAD=4nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=50_RCNRAD=4nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');

% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.25_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=0.5_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=1_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=2_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=3_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=5_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');
% simdata{length(simdata)+1}=load('LEEM_SimResults\Ein=80eV_Dose=10_RCNRAD=1nm_5nmX5nmEBeam_Univ=30nmX30nmX30nm_px=1nm.mat');

%% Read the data
clc;
% close all;

addpath('F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\ElectronExposures_SEM\LoweV_SEM\Analysis_Code');

datapath={};
datapath{length(datapath)+1}='F:\Documents and Settings\sbhattarai\My Documents\Research\LEE_Interactions\ElectronExposures_SEM\LoweV_SEM\BareResist\DataFiles\';

for i = 1:length(datapath)
    [fname_tmp,path_tmp]=uigetfile([datapath{i} '*.dat'],'multiselect','on');
end

if ~iscell(fname_tmp)
    fnames2{1}=fname_tmp;
else
    fnames2=fname_tmp;
end

%%% Now read the files
% clc;
% close all;

TVdata = {};
SLOW3data={};
SLOW5data={};

for i = 1:length(fnames2)
    fname=[path_tmp fnames2{i}];
    [dataID,outdata]=ReadData(fname);
    
    for j = 1:length(outdata) % get the TV data
        if strcmp(outdata{j}.sweepparm,'SR=TV2 (33 ms/frame)')==1
            idx=length(TVdata)+1;
            TVdata{idx}=outdata{j};
            TVdata{idx}.dataID=dataID;
        end
    end
    
    for j = 1:length(outdata) % get the SLOW3 data
        if strcmp(outdata{j}.sweepparm,'SR=SLOW3 (8 s/frame)')==1
            idx=length(SLOW3data)+1;
            SLOW3data{idx}=outdata{j};
            SLOW3data{idx}.dataID=dataID;
        end
    end
    
    for j = 1:length(outdata) % get the SLOW5 data
        if strcmp(outdata{j}.sweepparm,'SR=SLOW5 (30 s/frame)')==1
            idx=length(SLOW5data)+1;
            SLOW5data{idx}=outdata{j};
            SLOW5data{idx}.dataID=dataID;
        end
    end
end

pltstyle={'-ob','-ok','-or'};

data29={};
data29.dose=[];
data29.rt=[];

data49={};
data49.dose=[];
data49.rt=[];

data91={};
data91.dose=[];
data91.rt=[];

% fig1=figure;
% hold on;
% ylim([20 32]);
% set(gca,'XScale','log');

analysis_dataset=SLOW3data;

dosemin=Inf;
for i = 1:length(analysis_dataset)
    Iem=analysis_dataset{i}.contrcurve(:,1);
    texp=analysis_dataset{i}.contrcurve(:,2);
    Area=analysis_dataset{i}.contrcurve(:,3);
    rt=analysis_dataset{i}.contrcurve(:,4);
    dose_tmp=Iem.*texp./Area;
    dose_tmp=dose_tmp(rt~=0);
    rt=rt(rt~=0);
    
    switch analysis_dataset{i}.E
        case 29
            data29.dose=[data29.dose;dose_tmp];
            data29.rt=[data29.rt;rt];
            if min(data29.dose)<dosemin
                dosemin=min(data29.dose);
            end
%             plot(dose,rt,pltstyle{1},'linewidth',3.0);
        case 49
            data49.dose=[data49.dose;dose_tmp];
            data49.rt=[data49.rt;rt];
            if min(data49.dose)<dosemin
                dosemin=min(data49.dose);
            end
%             plot(dose,rt,pltstyle{2},'linewidth',3.0);
        case 91
            data91.dose=[data91.dose;dose_tmp];
            data91.rt=[data91.rt;rt];
            if min(data91.dose)<dosemin
                dosemin=min(data91.dose);
            end
%             plot(dose,rt,pltstyle{3},'linewidth',3.0);
    end
end

[data29.dose,idx]=sort(data29.dose);
data29.dose=data29.dose./dosemin;
data29.rt=data29.rt(idx);

[data49.dose,idx]=sort(data49.dose);
data49.dose=data49.dose./dosemin;
data49.rt=data49.rt(idx);

[data91.dose,idx]=sort(data91.dose);
data91.dose=data91.dose./dosemin;
data91.rt=data91.rt(idx);

%%% Just plot the raw data
% clc;
% close all;

% figure;
% plot(data29.dose,data29.rt,'-ob','linewidth',3.0,'MarkerSize',12);
% hold on;
% plot(data49.dose,data49.rt,'-ok','linewidth',3.0,'MarkerSize',12);
% plot(data91.dose,data91.rt,'-or','linewidth',3.0,'MarkerSize',12);
% xlabel('Relative Dose','fontsize',18);
% ylabel('Resist Thickness (nm)','fontsize',18);
% set(gca,'fontsize',20);

%% post-proc of LEEM analysis results: RCN/Diffn/Base Model [set up parameters]
clc;
% close all;

base_load=0.04; % PAG is 0.5/nm3;
acid_difflen=16; % nm, over a 90 s PEB time at 100 deg. C.
base_difflen=5; % nm, over a 90 s PEB time at 100 deg. C.
prot_load=1; % /nm3
kD=1.25; % Deprotection Rate (nm3/s)
kQ=5; % A/B Quenching Rate (nm3/s)
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
depr_thr=0.5;

%% post-proc of LEEM analysis results: RCN/Diffn/Base Model [Just get the rcn/diffn psf]
clc;
close all;

tmpdata=simdata{1}; % just a placeholder to get grids etc.
tmpdata.acidimg=zeros(tmpdata.univ.npx);
[nx,ny,nz]=size(tmpdata.acidimg);
xyz_tmp=round([nx,ny,nz]./2);
xo=xyz_tmp(1);yo=xyz_tmp(2);zo=xyz_tmp(3);
tmpdata.acidimg(xo,yo,zo)=2;

outparms=RcnDiffn(physparms,simparms,tmpdata);

deprimg=outparms.deprimg;
deprimg=deprimg./max(deprimg(:));

tmpidx=find(deprimg>=0.5);
[xidx,yidx,zidx]=ind2sub(size(deprimg),tmpidx);
radius=sqrt((xidx-xo).^2+(yidx-yo).^2+(zidx-zo).^2);

figure;hist(radius(:),50)

figure;imagesc(sum(deprimg,3));colorbar;

%% slice the deprotection image
clc;
% close all;

[x,y,z]=meshgrid(1:50,1:50,1:50);
depr_slice=slice(x,y,z,deprimg,[45 2 4],[52 52 50],36);
xlabel('x');
ylabel('y');
zlabel('z');

%% post-proc of LEEM analysis results: RCN/Diffn/Base Model [Run the actual sim]
% close all;

dose=[];
zvec2=[];

ntrials=3;
deprthr=0.5;

deprmean=[];deprmean2=[];deprmean3=[];
deprstd=[];deprstd2=[];deprstd3=[];
tstart=tic;
fig1=figure;
simparms.zadd_npx=0;
dose=zeros(1,length(simdata));
rcndiffndata={};
for i = 1:length(simdata)
% for i = [5 10]
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
        imagesc(mean(outparms.deprimg,3));caxis([0 1]);colorbar;
        drawnow;
    %     NdeprData=Count_thruz(outparms.deprimg,simdata{i}.elecimg_inc);
    %     zvec(i,1:length(simdata{i}.zvec))=simdata{i}.zvec;
        zvec=simdata{i}.univ.grid.z(:,1,1);
        zvec=zvec-zvec(1); % move such that top of resist is 0
    %     deprmean(i,1:length(NdeprData.mean))=NdeprData.mean;

    %     depth(i)=depth_calc(simdata{i}.zvec,NdeprData.mean,depr_thr);

    %     [depth(i),depthmat]=depth_calc2(outparms.deprimg,zvec,simdata{i}.elecimg_inc,0.5);
    %     [depth(i),depth2(i),depthmat]=depth_calc2(outparms.deprimg,zvec,simdata{i}.elecimg_inc,0.5);

    %     depr_thr=0.5;
        deprdata=depr_thruz(outparms,elecimg_bin,simdata{i}.acidimg,deprthr);

%         zvec2(i,1:length(zvec))=zvec';
%         deprmean(i,1:length(deprdata.mean))=deprdata.mean;
%         deprmean2(i,1:length(deprdata.mean))=deprdata.mean2;
%         deprstd(i,1:length(deprdata.std))=deprdata.mean;

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
    deprmean3(i,1:size(deprmean2_tmp,2))=mean(deprmean3_tmp,1);
    deprstd(i,1:size(deprmean_tmp,2))=std(deprmean_tmp,0,1); % std(X,0,DIM)
    deprstd2(i,1:size(deprmean2_tmp,2))=std(deprmean2_tmp,0,1); % std(X,0,DIM)
    deprstd3(i,1:size(deprmean2_tmp,2))=std(deprmean3_tmp,0,1); % std(X,0,DIM)
end
fprintf('Rcn/Diffn Calculations complete\n\n');
tend=toc(tstart);
fprintf('...%.2f s = %.2f min\n',tend,tend/60);

% figure;
% plot(dose,50-depth,'-o')
% set(gca,'XScale','log');
% xlabel('Dose');
% ylabel('Depth');
% xlim([0.1 100]);
% ylim([22 32])

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

%% Extract thickness loss at the threshold from the deprmean matrix above
% clc;
% close all;

rtlost=[];
rtlost_std=[];
depr_thr=0.3;
analysis_entity=deprmean2;

count=1;
for i = 1:size(analysis_entity,1)
% for i = [1 2 3 5 6 7 8 9 11 12 13 14 15 17 18]
    ztmp=zvec2(i,:);
    ztmp=[-simparms.zadd_npx:1:ztmp(1)-1 ztmp];
%     ztmp=ztmp-simparms.zadd_npx; % nmpp=1 by default, double check in future
    deprtmp=analysis_entity(i,:);
    deprtmp_std=deprstd(i,:);
    
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
    end
    rtlost(count)=mean(rt_tmp);
    rtlost_std(count)=std(rt_tmp);
    count=count+1;
end

% clc;

t_unexp=50;
rt_ref=24;
rt=t_unexp-rtlost;

%%%%%% rt vs. logarithmic dose fit
% p=polyfit(log10(dose),rt,1);
% p

% figure;plot(dose(1:4),50-rtlost(1:4),'-ob','linewidth',3.0);
% hold on;
% plot(dose(5:8),50-rtlost(5:8),'-ok','linewidth',3.0);
% plot(dose(9:12),50-rtlost(9:12),'-or','linewidth',3.0)
% set(gca,'XScale','log');

% figure;plot(dose,rt,'-ob');
% set(gca,'XScale','log');
% figure;errorbar(dose,t_unexp-rtlost,rtlost_std,'-ob')

%% create the data structures that hold dose, thickness values
sim29={};
sim49={};
sim91={};
sim29.dose=dose(1,1:5);
sim29.rt=t_unexp-rtlost(1,1:5);
sim29.rt_eb=rtlost_std(1,1:5);
% try
%     sim49.dose=dose(1,6:10);
%     sim49.rt=t_unexp-rtlost(1,6:10);
%     sim49.rt_eb=rtlost_std(1,6:10);
% catch
% %     fprintf('error in sim49');
% end

try
    sim91.dose=dose(1,6:10);
    sim91.rt=t_unexp-rtlost(1,6:10);
    sim91.rt_eb=rtlost_std(1,6:10);
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

p1=polyfit(log10(sim29.dose(3:5)),sim29.rt(3:5),1);
p2=polyfit(log10(sim91.dose(3:5)),sim91.rt(3:5),1);

figure;
set(gcf,'Position',gcf_pos);
set(gca,'Position',gca_pos);
errorbar(1.*sim29.dose,sim29.rt,sim29.rt_eb,'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
% plot(sim29.dose,polyval(p1,log10(sim29.dose)),'k','linewidth',3.0);
% errorbar(1.*sim49.dose,sim49.rt,sim49.rt_eb,'-ok','linewidth',3.0,'MarkerSize',12);
errorbar(sim91.dose,sim91.rt,sim91.rt_eb,'-or','linewidth',3.0,'MarkerSize',12);
% plot(sim91.dose,polyval(p2,log10(sim91.dose)),'k','linewidth',3.0);
xlabel('Dose (e^-/nm^2)','fontsize',30);
ylabel('Thicknss (nm)','fontsize',30);
set(gca,'fontsize',30,'linewidth',3.0,'XScale','log');
xlim([0.5 20]);
ylim([38 50]);
legend('29 eV','49 eV','91 eV');
legend('29 eV','91 eV');

% figure;
% set(gcf,'Position',gcf_pos);
% set(gca,'Position',gca_pos);
% errorbar(0.625.*sim29.dose,sim29.rt,sim29.rt_eb,'-ob','linewidth',3.0,'MarkerSize',12);
% hold on;
% % errorbar(49/91.*sim49.dose,sim49.rt,sim49.rt_eb,'-ok','linewidth',3.0,'MarkerSize',12);
% errorbar(sim91.dose,sim91.rt,sim91.rt_eb,'-or','linewidth',3.0,'MarkerSize',12);
% xlabel('Dose (e^-/nm^2)','fontsize',30);
% ylabel('Thicknss (nm)','fontsize',30);
% set(gca,'fontsize',30,'linewidth',3.0,'XScale','linear');
% % xlim([0.7 300]);
% % ylim([40 52]);
% legend('29 eV','49 eV','91 eV');

%% plot means from other trials
clc;
close all;

t1=load('data_Dose=[1,2,4,8]_E=[29,49,91]_RcnDiffn_T1.mat');
t2=load('data_Dose=[1,2,4,8]_E=[29,49,91]_RcnDiffn_T2.mat');
t3=load('data_Dose=[1,2,4,8]_E=[29,49,91]_RcnDiffn_T3.mat');

rt1=[t1.sim29.rt;t2.sim29.rt;t3.sim29.rt];
rt2=[t1.sim49.rt;t2.sim49.rt;t3.sim49.rt];
rt3=[t1.sim91.rt;t2.sim91.rt;t3.sim91.rt];

figure;
plot([1 2 4 8],mean(rt1,1),'-ob');
hold on;
plot([1 2 4 8],mean(rt2,1),'-ok');
plot([1 2 4 8],mean(rt3,1),'-or');
set(gca,'XScale','linear');

%% # of acids in the LEEM simulations
clc;
close all;

nAcids=[];
for i = 1:length(simdata)
    nAcids(i)=sum(simdata{i}.acidimg(:));
end

figure;
plot(nAcids,'-ob');

%% 29 eV vs. 91 eV 
% clc;
close all;

figure;
% plot(deprmean2(1,:),'-ob','linewidth',3.0);hold on;
errorbar(deprmean2(1,:),deprstd2(1,:),'-ob','linewidth',2.0);hold on;
% plot(deprmean2(3,:),'-ok','linewidth',3.0);
errorbar(deprmean2(3,:),deprstd2(3,:),'-ok','linewidth',2.0);
xlabel('z (nm)');
ylabel('Deprotection Level');
legend('29 eV (4 e^-/nm^2)','91 eV (4 e^-/nm^2');
set(gca,'fontsize',30,'linewidth',3.0);
ylim([0 0.8]);
xlim([0 50]);

figure;
% plot(deprstd(1,:),'-ob','linewidth',3.0);hold on;
errorbar(deprmean2(2,:),deprstd2(2,:),'-ob','linewidth',2.0);hold on;
% plot(deprstd(3,:),'-ok','linewidth',3.0);
errorbar(deprmean2(3,:),deprstd2(3,:),'-ok','linewidth',2.0);
xlabel('z (nm)');
ylabel('Deprotection Level');
legend('29 eV (12 e^-/nm^2)','91 eV (4 e^-/nm^2');
set(gca,'fontsize',30,'linewidth',3.0);
ylim([0 0.8]);
xlim([0 50]);

%% Plot thru-z deprotection levels for publication purposese

clc;
close all;
gcf_pos=[680   301   851   677];
gca_pos=[0.1300    0.1100    0.7750    0.8150];

figure;
plot(deprmean2([1 :5],:)','linewidth',3.0);
title('29 eV');
xlim([0 20]);ylim([0 1]);
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);

figure;
plot(deprmean2([6 :10],:)','linewidth',3.0);
title('49 eV');
xlim([0 20]);ylim([0 1]);
set(gca,'fontsize',30,'linewidth',3.0);
set(gca,'Position',gca_pos);
set(gcf,'Position',gcf_pos);

% figure;
% plot(deprmean2([13 :18],:)','linewidth',3.0);
% title('91 eV');
% xlim([0 20]);ylim([0 1]);
% set(gca,'fontsize',30,'linewidth',3.0);
% set(gca,'Position',gca_pos);
% set(gcf,'Position',gcf_pos);

%% Plotting acid data for publication purposes

% close all;
selidx=3;
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
outparms=rcndiffndata{5};
tmp=outparms.deprimg(:,:,23:29);
% tmp(outparms.protimg_pre(:,:,23:29)==0)=1;
tmp2=mean(tmp,3);
tmp3=tmp2;
% tmp3=[ones(1,size(tmp3,2));tmp3];
tmp3(tmp3>=thr_tmp)=1;
tmp3(tmp3<thr_tmp)=0;

figure;
imagesc(tmp3);colorbar;

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

figure;
imagesc(tmp3);colorbar;
% imagesc(edge(tmp3));colorbar;

%% plot thru-z acid and deprotection image cross-sections for publications 
clc;
% close all;
gcf_pos=[370   320   870   658];
gca_pos=[0.1300    0.1337    0.7750    0.7913];
acidsum=[];
selidx=[7,8];
selidx=[1,2];
    acidsum=[];

for j = 1:length(selidx)
    acidtmp=simdata{selidx(j)}.acidimg;
    ztmp=simdata{selidx(j)}.univ.grid.z;
    ztmp=ztmp(:,1,1);
    ztmp=ztmp(:);

    for i = 1:size(acidtmp,1)
        tmpmat=acidtmp(i,:,:);
        acidsum(i,j)=sum(tmpmat(:));
    end
end

figure(fig1);
% figure
hold on;
plot(ztmp,mean(acidsum,2),'-.b','linewidth',3.0,'MarkerSize',12)

%% publication plots: Total # of acids vs. dose vs. energy
% clc;
close all;

gcf_pos=[370   320   870   658];
gca_pos=[0.1300    0.1337    0.7750    0.7913];

nacids=[];
dosetmp=[];
acidsum=[];
for i = 1:length(simdata)
    acidimg=simdata{i}.acidimg;
%     acidimg=acidimg(1:5,:,:);
    dosetmp(i,1)=simdata{i}.Dose;
    for j = 1:size(acidimg,1)
        acidtmp=acidimg(j,:,:);
        nacids(i,j)=sum(acidtmp(:));
    end
    acidimg=acidimg(1:5,:,:);
    acidsum(i)=sum(acidimg(:));
end

acids29=nacids(1:8,:);
dose29=dosetmp(1:8);

acids49=nacids(9:14,:);
dose49=dosetmp(9:14);
[dose49,idx]=sort(dose49);
acids49=acids49(idx,:);

acids91=nacids(15:end,:);
dose91=dosetmp(15:end);

figure;
plot([6 25 100],[acidsum(1) acidsum(5) acidsum(8)],'-ob','linewidth',3.0); % 29 eV
hold on;
plot([6 25 100],[acidsum(9) acidsum(11) acidsum(13)],'-ok','linewidth',3.0); % 49 eV
plot([6 25 100],[acidsum(15) acidsum(17) acidsum(19)],'-or','linewidth',3.0);


figure;
plot([6 25 100],[acidsum(1) acidsum(5) acidsum(8)]./acidsum(19),'-ob','linewidth',3.0); % 29 eV
hold on;
plot([6 25 100],[acidsum(9) acidsum(11) acidsum(13)]./acidsum(19),'-ok','linewidth',3.0); % 49 eV
plot([6 25 100],[acidsum(15) acidsum(17) acidsum(19)]./acidsum(19),'-or','linewidth',3.0);

%%% plot thru-z acids using numbers from above cell
% figure;
% hold on;
% plot(dose29,acids29,'--ob','linewidth',3.0,'MarkerSize',12);
% hold on;
% plot(dose49,acids49,'--ok','linewidth',3.0,'MarkerSize',12);
% plot(dose91,acids91,'--or','linewidth',3.0,'MarkerSize',12);
% xlabel('Dose');
% ylabel('# of Acids');
% legend('29 eV','49 eV','91 eV');

%%% dose = 6;
figure;
plot(zvec',nacids(1,:),'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
plot(zvec',nacids(9,:),'-ok','linewidth',3.0,'MarkerSize',12);
plot(zvec',nacids(15,:),'-or','linewidth',3.0,'MarkerSize',12);
xlim([0 15]);ylim([0 50]);
xlabel('Z (nm)');
ylabel('No. of Acids');
set(gca,'fontsize',30,'linewidth',3.0);

%%% dose = 25;
figure;
plot(zvec',nacids(5,:),'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
plot(zvec',nacids(11,:),'-ok','linewidth',3.0,'MarkerSize',12);
plot(zvec',nacids(17,:),'-or','linewidth',3.0,'MarkerSize',12);
xlim([0 15]);ylim([0 50]);
xlabel('Z (nm)');
ylabel('No. of Acids');
set(gca,'fontsize',30,'linewidth',3.0);

%%% dose = 100;
figure;
plot(zvec',nacids(8,:),'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
plot(zvec',nacids(13,:),'-ok','linewidth',3.0,'MarkerSize',12);
plot(zvec',nacids(19,:),'-or','linewidth',3.0,'MarkerSize',12);
xlim([0 15]);ylim([0 50]);
xlabel('Z (nm)');
ylabel('No. of Acids');
set(gca,'fontsize',30,'linewidth',3.0);


% figure;
% plot(nacids(1,:)./nacids(1,2),'-ob','linewidth',3.0,'MarkerSize',12);
% hold on;
% plot(nacids(9,:)./nacids(9,2),'-ok','linewidth',3.0,'MarkerSize',12);
% plot(nacids(15,:)./nacids(15,2),'-or','linewidth',3.0,'MarkerSize',12);
% xlim([0 15]);
% 
% figure;
% plot(nacids(5,:)./nacids(5,3),'-ob','linewidth',3.0,'MarkerSize',12);
% hold on;
% plot(nacids(11,:)./nacids(11,3),'-ok','linewidth',3.0,'MarkerSize',12);
% plot(nacids(17,:)./nacids(17,3),'-or','linewidth',3.0,'MarkerSize',12);
% xlim([0 15]);
% 
% figure;
% plot(nacids(8,:)./nacids(8,3),'-ob','linewidth',3.0,'MarkerSize',12);
% hold on;
% plot(nacids(13,:)./nacids(13,3),'-ok','linewidth',3.0,'MarkerSize',12);
% plot(nacids(19,:)./nacids(19,3),'-or','linewidth',3.0,'MarkerSize',12);
% xlim([0 15]);


%%% 29 eV dose=6
idx=5;
acidimg=simdata{idx}.acidimg;
acidimg=acidimg(1:5,:,:);
density29=sum(acidimg(:))/(9*9*5); % 9X9nm2 determined by looking at images
figure;imagesc(reshape(sum(acidimg,1),50,50));colorbar;
title(['29 eV; dose = ' num2str(dose(idx))]);

%%% 49 eV dose=6
idx=11;
acidimg=simdata{idx}.acidimg;
acidimg=acidimg(1:5,:,:);
density49=sum(acidimg(:))/(9*9*5); % 9X9nm2 determined by looking at images
figure;imagesc(reshape(sum(acidimg,1),50,50));colorbar;
title(['49 eV; dose = ' num2str(dose(idx))]);

%%% 91 eV dose=6
idx=17;
acidimg=simdata{idx}.acidimg;
acidimg=acidimg(1:5,:,:);
density91=sum(acidimg(:))/(9*9*5); % 9X9nm2 determined by looking at images
figure;imagesc(reshape(sum(acidimg,1),50,50));colorbar;
title(['91 eV; dose = ' num2str(dose(idx))]);

[density29,density49,density91]

%% plot thru-z deprotection profiles
% clc;
close all;

ztmp=simdata{1}.univ.grid.z(:,1,1);
 
gcf_pos=[680   223   958   755];
gca_pos=[0.1300    0.1166    0.7750    0.8084];

figure;
% errorbar(ztmp',deprmean(11:end),deprstd(11:end),'-ob')
% errorbar(ztmp-min(ztmp),deprmean2(1,:)',deprstd2(1,:)','-ob','linewidth',3.0,'MarkerSize',6);
plot(ztmp-min(ztmp),deprmean2(6,:)','-ob','linewidth',3.0,'MarkerSize',12);
hold on;
% errorbar(ztmp-min(ztmp),deprmean2(2,:)',deprstd2(2,:)','-ok','linewidth',3.0,'MarkerSize',6);
plot(ztmp-min(ztmp),deprmean2(7,:)','-ok','linewidth',3.0,'MarkerSize',12);
% errorbar(ztmp-min(ztmp),deprmean2(3,:)',deprstd2(2,:)','-or','linewidth',3.0,'MarkerSize',6);
plot(ztmp-min(ztmp),deprmean2(8,:)','-or','linewidth',3.0,'MarkerSize',12);
% errorbar(ztmp-min(ztmp),deprmean2(4,:)',deprstd2(2,:)','-og','linewidth',3.0,'MarkerSize',6);
plot(ztmp-min(ztmp),deprmean2(9,:)','-og','linewidth',3.0,'MarkerSize',12);
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

%% append other data to existing data structures from above
% clc;
close all;

tmpdata=load('Data_Dose=200_10Trials.mat');

sim29_2.dose=[sim29.dose tmpdata.sim29.dose];
sim29_2.rt=[sim29.rt tmpdata.sim29.rt];
sim29_2.rt_eb=[sim29.rt_eb tmpdata.sim29.rt_eb];

sim49_2.dose=[sim49.dose tmpdata.sim49.dose];
sim49_2.rt=[sim49.rt tmpdata.sim49.rt];
sim49_2.rt_eb=[sim49.rt_eb tmpdata.sim49.rt_eb];

sim91_2.dose=[sim91.dose tmpdata.sim91.dose];
sim91_2.rt=[sim91.rt tmpdata.sim91.rt];
sim91_2.rt_eb=[sim91.rt_eb tmpdata.sim91.rt_eb];

%% add stats if you have other simulation runs data
clc;
close all;

% load('ResSim_data\dataset2.mat');

adddata=load('ResSim_data\data_29eV_T2_thr=0.4.mat');

tmpdata.dose=adddata.dose;
tmpdata.rt=adddata.t_unexp-adddata.rtlost;

sim29.rt=sim29.rt+adddata.t_unexp-t_unexp;

tmpdata2=AddStat(sim29,tmpdata);
sim29=tmpdata2;

%% optimizations and comparing with real data:

fun=@(dose_ref,dose,rt,rt_ref) abs(interp1(dose,rt,dose_ref,'linear')-rt_ref);

[minval,idx]=min(abs(sim29.rt-rt_ref));dose_ref_sim=mean(sim29.dose(idx));
dose_ref_sim=mean(sim29.dose);
dose_ref_sim=fminsearch(@(dose_ref_sim) fun(dose_ref_sim,sim29.dose,sim29.rt,rt_ref),dose_ref_sim);
sim29.doseref=dose_ref_sim;

[minval,idx]=min(abs(sim49.rt-rt_ref));dose_ref_sim=mean(sim49.dose(idx));
dose_ref_sim=mean(sim49.dose);
dose_ref_sim=fminsearch(@(dose_ref_sim) fun(dose_ref_sim,sim49.dose,sim49.rt,rt_ref),dose_ref_sim);
sim49.doseref=dose_ref_sim;

[minval,idx]=min(abs(sim91.rt-rt_ref));dose_ref_sim=mean(sim91.dose(idx));
% dose_ref_sim=median(sim91.dose);
dose_ref_sim=fminsearch(@(dose_ref_sim) fun(dose_ref_sim,sim91.dose,sim91.rt,rt_ref),dose_ref_sim);
idx1=find(sim91.rt>=rt_ref);if ~isempty(idx1) idx1=idx1(end);end
idx2=find(sim91.rt<rt_ref);if ~isempty(idx2) idx2=idx2(1);end
p=polyfit(sim91.dose(idx1:idx2),sim91.rt(idx1:idx2),1);dose_ref_sim=(rt_ref-p(2))/p(1);
sim91.doseref=dose_ref_sim;

[minval,idx]=min(abs(data29.rt2_mean-rt_ref));dose_ref_expt=mean(data29.dose(idx));
dose_ref_expt=fminsearch(@(dose_ref_expt) fun(dose_ref_expt,data29.dose2,data29.rt2_mean,rt_ref),dose_ref_expt);
data29.doseref=dose_ref_expt;

[minval,idx]=min(abs(data49.rt2_mean-rt_ref));dose_ref_expt=mean(data49.dose(idx));
dose_ref_expt=fminsearch(@(dose_ref_expt) fun(dose_ref_expt,data49.dose2,data49.rt2_mean,rt_ref),dose_ref_expt);
data49.doseref=dose_ref_expt;

[minval,idx]=min(abs(data91.rt2_mean-rt_ref));dose_ref_expt=mean(data91.dose(idx));
dose_ref_expt=fminsearch(@(dose_ref_expt) fun(dose_ref_expt,data91.dose2,data91.rt2_mean,rt_ref),dose_ref_expt);
data91.doseref=dose_ref_expt;

figure;
plot(sim29.dose./sim29.doseref,sim29.rt,'-b','linewidth',3.0);
errorbar(sim29.dose./sim29.doseref,sim29.rt,sim29.rt_eb,'-b','linewidth',3.0);
hold on;
plot(data29.dose2./data29.doseref,data29.rt2_mean,'ok','linewidth',3.0,'MarkerSize',12);
set(gca,'XScale','log');
title('29 eV');
xlabel('Dose');
ylabel('rt (nm)');
set(gca,'XScale','log','YScale','linear');xlim([0.05 100]);ylim([12 32])

figure;
plot(sim49.dose./sim49.doseref,sim49.rt,'-b','linewidth',3.0);
hold on;
plot(data49.dose2./data49.doseref,data49.rt2_mean,'ok','linewidth',3.0,'MarkerSize',12);
set(gca,'XScale','log');
title('49 eV');
xlabel('Dose');
ylabel('rt (nm)');
set(gca,'XScale','log','YScale','linear');xlim([0.05 100]);ylim([12 32])

figure;
plot(sim91.dose./sim91.doseref,sim91.rt,'--b','linewidth',3.0);
hold on;
plot(data91.dose2./data91.doseref,data91.rt2_mean,'or','linewidth',3.0,'MarkerSize',12);
set(gca,'XScale','log');
title('91 eV');
xlabel('Dose');
ylabel('rt (nm)');
set(gca,'XScale','log','YScale','linear');xlim([0.05 100]);ylim([12 32])

figure;
plot(sim29.dose,sim29.rt,'b','linewidth',3.0);
hold on;
plot(sim49.dose,sim49.rt,'k','linewidth',3.0);
plot(sim91.dose,sim91.rt,'r','linewidth',3.0);
xlabel('Dose (electrons/pixel');
ylabel('Thickenss (nm)');
legend('29 eV','49 eV','91 eV');
set(gca,'XScale','log','YScale','linear');xlim([0.05 100]);ylim([12 32])

%% Single convolution model [No 3-D SuMMIT model]
clc;
close all;

Nacids=[];
Nacids_per_px=[];
for i = 1:length(simdata)
    acidimg=simdata{i}.acidimg;
    outdata=getNthruz(acidimg);
    Nacids(:,i)=outdata.N;
    Nacids_per_px(:,i)=outdata.N_per_px;
end

idx_29=[1 8];

gauss_sigma=8;
flt1D=fspecial('gaussian',[size(Nacids,1)*3 1],gauss_sigma);
Nacids_filt=[];
count=1;
for i = idx_29(1):idx_29(2)
    vec=Nacids(:,i);
    vecf=conv(vec,flt1D,'same');
    Nacids_filt(:,count)=vecf;
    count=count+1;
end

%% post-proc of LEEM results [dose sweep, only plot sim results]
clc;
close all;

dose1=[1 2 3 5 10 20]; % 15 eV
dose2=[0.25 0.5 1 2 3 5 10 20]; % 20 eV
dose3=[0.25 0.5 1 2 3 5 10 20 50]; % 40 eV
dose4=[0.1 0.25 0.5 1 2 3 5 10 20 50]; % 80 eV

dose_sim=[];
depth_sim=[];

dose_sim(size(dose_sim,1)+1,1:length(dose1))=dose1;
dose_sim(size(dose_sim,1)+1,1:length(dose2))=dose2;
dose_sim(size(dose_sim,1)+1,1:length(dose3))=dose3;
dose_sim(size(dose_sim,1)+1,1:length(dose4))=dose4;

depth_sim(size(depth_sim,1)+1,1:length(dose1))=depth(size(depth_sim,2)+1:size(depth_sim,2)+length(dose1));
depth_sim(size(depth_sim,1)+1,1:length(dose2))=depth(1+length(dose1):length(dose1)+length(dose2));
depth_sim(size(depth_sim,1)+1,1:length(dose3))=depth(1+length(dose1)+length(dose2):length(dose1)+length(dose2)+length(dose3));
depth_sim(size(depth_sim,1)+1,1:length(dose4))=depth(1+length(dose1)+length(dose2)+length(dose3):length(dose1)+length(dose2)+length(dose3)+length(dose4));


tmpidx=find(dose_sim==0);
dose_sim(tmpidx)=NaN;
depth_sim(tmpidx)=NaN;

%%%% Dose-Normalized Plot
figure;
plot(15/80.*dose_sim(1,:),depth_sim(1,:),'-ob','linewidth',3.0);
hold on;
plot(20/80.*dose_sim(2,:),depth_sim(2,:),'-or','linewidth',3.0);
plot(40/80.*dose_sim(3,:),depth_sim(3,:),'-ok','linewidth',3.0);
plot(dose_sim(4,:),depth_sim(4,:),'-og','linewidth',3.0);
xlabel('Energy-Normalized Dose','fontsize',30);
ylabel('Thickness Loss (nm)','fontsize',30);
set(gca,'XScale','log','fontsize',30,'linewidth',3.0);
box on;

%%%% Absolute-Value Plot
figure;
plot(dose_sim(1,:),depth_sim(1,:),'-ob','linewidth',3.0);
hold on;
plot(dose_sim(2,:),depth_sim(2,:),'-or','linewidth',3.0);
plot(dose_sim(3,:),depth_sim(3,:),'-ok','linewidth',3.0);
plot(dose_sim(4,:),depth_sim(4,:),'-og','linewidth',3.0);
xlabel('Dose (e^-/nm^2)','fontsize',30);
ylabel('Thickness Loss (nm)','fontsize',30);
legend('15 eV','20 eV','40 eV','80 eV');
set(gca,'XScale','log','fontsize',30,'linewidth',3.0);
box on;

%% post-proc of LEEM analysis results [when swept parameter is dose]
% clc;
close all;

%%%% pick out the simulation data
Esim=[40 80];
% Esim=20;
% E_expt=expt_data.E;
E_expt=[40 60 80];
% E_expt=[15 20 25 30];

% dose1=[0.5 1 2 3 5 10 20 50];
% dose1=[0.25 0.5 1 2 3 5 10 20];
dose1=[0.25 0.5 1 2 3 5 10 20 50];
dose2=[0.1 0.25 0.5 1 2 3 5 10 20 50];
% dose2=[0.1 0.25 0.5 1 2 3];
% dose2=[];
dose_sim=[];
depth_sim=[];

dose_sim(size(dose_sim,1)+1,1:length(dose1))=dose1;
depth_sim(size(depth_sim,1)+1,1:length(dose1))=depth(size(depth_sim,2)+1:size(depth_sim,2)+length(dose1));

if ~isempty(dose2)
    dose_sim(size(dose_sim,1)+1,1:length(dose2))=dose2;
    depth_sim(size(depth_sim,1)+1,1:length(dose2))=depth(size(depth_sim,2)+1:size(depth_sim,2)+length(dose2));
end

idx=find(dose_sim==0);
dose_sim(idx)=NaN;
depth_sim(idx)=NaN;

%%%% load the experimental data tables
expt_data_path='ExptData\';
expt_data_fname='LEEMSummary_E=[40,60,80].mat';
% expt_data_fname='LEEMSummary_E=[15,20,25,30].mat';
expt_data=load([expt_data_path expt_data_fname]);

flag=ismember(E_expt,Esim);
idx=find(flag==1);

fig1=figure;hold on;
% figure(fig1);plot(dose_sim',depth_sim','-','MarkerSize',12);

simdata_count=1;
sim_plt_styles={'-b','-r','-k','-g','--b','--r','--k','--g'};
sim_plt_styles=repmat(sim_plt_styles,[1,100]);
% sim_plt_styles=sim_plt_styles(2:end);

expt_plt_styles={'ob','or','ok','og'};
expt_plt_styles=repmat(expt_plt_styles,[1,100]);
% expt_plt_styles=expt_plt_styles(2:end);
    
for idx2=idx
% for idx2=2
    %%%% Experiment vectors
%     if idx2==1 % just a simple hack, make more elegant in future
        dose_expt=expt_data.data{idx2}.reldose;
        rt_expt=expt_data.rt;
        tloss_expt=rt_expt-expt_data.data{idx2}.thickness;
%     end
    
    %%%% Simulation vectors
    dose_sim2=dose_sim(simdata_count,:);
    tloss_sim=depth_sim(simdata_count,:);
    
    %%%% find the anchoring point in expt data
    anchoridx=find(tloss_expt>min(tloss_sim));
    if idx2==1
        dose_anchor=dose_expt(anchoridx(2));
        tloss_anchor=tloss_expt(anchoridx(2));
    else
        dose_anchor=dose_expt(anchoridx(2));
        tloss_anchor=tloss_expt(anchoridx(2));
    end
    
    %%%% interpolate and find the corresponding point in sim data
    anchoridx=find(tloss_sim==tloss_anchor);
    if ~isempty(anchoridx)
        dose_sim2=dose_sim2./dose_sim2(anchoridx);
    else
        anchoridx_2=find(tloss_sim<tloss_anchor);
        if ~isempty(anchoridx_2)
            anchoridx_2=anchoridx_2(end);
            anchoridx_3=find(tloss_sim>tloss_anchor);
            anchoridx_3=anchoridx_3(1);
            p=polyfit(dose_sim2(anchoridx_2:anchoridx_3),tloss_sim(anchoridx_2:anchoridx_3),1);
            dose_anchor_sim=(tloss_anchor-p(2))/p(1);
            
        else
            anchoridx_3=find(tloss_sim>tloss_anchor);
            anchoridx_3=anchoridx_3(1);
            p=polyfit((dose_sim2(anchoridx_3:anchoridx_3+1)),tloss_sim(anchoridx_3:anchoridx_3+1),1);
            dose_anchor_sim=(tloss_anchor-p(2))/p(1);
%             dose_anchor_sim=10^dose_anchor_sim;
        end
        dose_sim2=dose_sim2./dose_anchor_sim.*dose_anchor;
%         dose_sim2=dose_sim2+dose_anchor-dose_anchor_sim;
    end
    
    dose_interp=dose_sim2(~isnan(dose_sim2));
    tloss_interp=tloss_sim(~isnan(dose_sim2));
    error=tloss_expt-interp1(dose_interp,tloss_interp,dose_expt,'nearest','extrap');
    rmsd=sqrt(mean(error.^2));
    fprintf('RMSD @ %d eV = %.4f\n',expt_data.data{idx2}.E,rmsd);
    
    figure(fig1);
%     plot(dose_sim(simdata_count,:),depth_sim(simdata_count,:),sim_plt_styles{simdata_count},'linewidth',3.0);
    plot(dose_sim2,tloss_sim,sim_plt_styles{simdata_count},'linewidth',3.0);
%     if idx2==1
        plot(dose_expt,tloss_expt,expt_plt_styles{simdata_count},'MarkerSize',20,'linewidth',3.0);
%     end
    xlabel('Rel. Dose','fontsize',30);
    ylabel('Thickness Lost (nm)','fontsize',30);
    set(gca,'XScale','log','fontsize',30,'linewidth',3.0);box on;
    ylim([0 16])
%     legend('40 eV Sim','40 eV expt','80 eV Sim','80 eV expt');
%     legend('40 eV Sim (rcnrad=1)','40 eV expt','40 eV Sim (rcnrad=2)');
    simdata_count=simdata_count+1;
end
fprintf('\n');
legend('40 eV Sim','40 eV Expt','80 eV Sim','80 eV Expt','Orientation','Horizontal');

% figure;%plot(dose_sim',depth_sim','-o','linewidth',3.0);set(gca,'XScale','log')
% plot(0.5.*dose_sim(1,:),depth_sim(1,:),'-ob','linewidth',3.0,'MarkerSize',12);
% % hold on; plot(dose_sim(2,:),depth_sim(2,:),'-or','linewidth',3.0,'MarkerSize',12);
% set(gca,'XScale','log')

%% TESTS: Analysis on slices thru-Z [Mean deprot. level, RMSD]

clc;close all;

deprimg=outparms.deprimg;

rmsd=[];slice_mean=[];
thr_min=0.1;
fig1=figure;
for i = 1:size(deprimg,1)
    
    slice=deprimg(i,:,:);
    slice=reshape(slice,[size(slice,2) size(slice,3)]);
    slice(slice<thr_min)=0;
    
    %%%% plot the image
    figure(fig1);
    imagesc(slice);colorbar;caxis([0 1]);
    title(sprintf('Slice %d',i));
    drawnow;
    
    rmsd(i)=sqrt(mean((slice(:)-0.5).^2));
    slice_mean(i)=mean(slice(slice>0));
end

figure;
plot(rmsd,'-ob','linewidth',3.0);
xlabel('Slice #');
ylabel('RMSD');

figure;
plot(slice_mean,'-ob','linewidth',3.0);
xlabel('Slice #');
ylabel('Mean');

%% TESTS: Analysis on slices thru-Z [1-D depth map]

clc;close all;

deprimg=outparms.deprimg;
depr_thr=0.5;
for i = 1:size(deprimg,2)
    for j = size(deprimg,3):-1:1
        deprvec=deprimg(:,i,j);
        deprvec(deprvec>=depr_thr)=1;
        deprvec(deprvec<depr_thr)=0;
        idx=find(deprvec==1);
        if ~isempty(idx)
            depth(i,j)=idx(end);
        else
            depth(i,j)=0;
        end
        
    end
end




