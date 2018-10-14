%% Read saved LEEM simulation files
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

%% Plotting LEEM data for publication purposes

close all;
selidx=8;
figure;
imagesc(sum(simdata{selidx}.acidimg,3));
colorbar;
xlim([10 40]);
ylim([1 25]);
% caxis([0 10]);
set(gca,'fontsize',30,'linewidth',3.0);
simdata{selidx}.Dose

%% load the data file
clear;
clc;
close all;

% load('ResSim_data\dataset2.mat');
load('ResSim_data\SavedDataSet_10Trials_04172017.mat');


%% plot thru-z acid and deprotection image cross-sections for publications 
% clc;
% close all;
gcf_pos=[370   320   870   658];
gca_pos=[0.1300    0.1337    0.7750    0.7913];
acidsum=[];
selidx=[14];
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
plot(ztmp,mean(acidsum,2),'-.k','linewidth',3.0,'MarkerSize',12)

%% plot the deprotection profiles for publications 
clc;
close all;

tmp=outparms.deprimg(:,:,23:29);

figure;
imagesc(mean(tmp,3));
colorbar;
xlim([10 35]);
ylim([1 25]);
caxis([0 1]);
set(gca,'fontsize',30,'linewidth',3.0);

%% Extract thickness loss at the threshold from the deprmean matrix above
clc;
% close all;

rtlost=[];
rtlost_std=[];
depr_thr=0.5;
analysis_entity=deprmean;
for i = 1:size(analysis_entity,1)
% for i = 16
    ztmp=zvec2(i,:);
    ztmp=[-simparms.zadd_npx:1:ztmp(1)-1 ztmp];
%     ztmp=ztmp-simparms.zadd_npx; % nmpp=1 by default, double check in future
    deprtmp=analysis_entity(i,:);
    deprtmp_std=deprstd(i,:);
    
    ntrials2=100;
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
    rtlost(i)=mean(rt_tmp);
    rtlost_std(i)=std(rt_tmp);
end

%%%
clc;

t_unexp=50;
rt_ref=24;

% figure;plot(dose,t_unexp-rtlost,'-ob');
% figure;errorbar(dose,t_unexp-rtlost,rtlost_std,'-ob')

sim29={};
sim49={};
sim91={};
% sim29.dose=dose(1,1:13);
% sim29.rt=t_unexp-rtlost(1,1:13);
% sim29.rt_eb=rtlost_std(1,1:13);
[sim29.dose,sim29.rt,sim29.rt_eb]=getMean(dose(1,1:13),t_unexp-rtlost(1,1:13),rtlost_std(1,1:13));

try
%     sim49.dose=dose(1,14:18);
%     sim49.rt=t_unexp-rtlost(1,14:18);
%     sim49.rt_eb=rtlost_std(1,14:18);
    [sim49.dose,sim49.rt,sim49.rt_eb]=getMean(dose(1,14:18),t_unexp-rtlost(1,14:18),rtlost_std(1,14:18));
catch
%     fprintf('error in sim49');
end

try
%     sim91.dose=dose(1,19:end);
%     sim91.rt=t_unexp-rtlost(1,19:end);
%     sim91.rt_eb=rtlost_std(1,19:end);
    [sim91.dose,sim91.rt,sim91.rt_eb]=getMean(dose(1,19:end),t_unexp-rtlost(1,19:end),rtlost_std(1,19:end));
catch
end

%%% find where rtlost==tunexposed
% simulation:
rt_sim=t_unexp-rtlost(1:end);
dose_sim=dose(1:end);

%% Analyze only the simulation results
clc;
% close all;

gcf_pos=[680   307   824   671];
gca_pos=[0.1300    0.1759    0.7750    0.7491];

figure;
errorbar(0.27.*sim29.dose,sim29.rt,sim29.rt_eb,'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
errorbar(0.4.*sim49.dose,sim49.rt,sim49.rt_eb,'-ok','linewidth',3.0,'MarkerSize',12);
errorbar(sim91.dose,sim91.rt,sim91.rt_eb,'-or','linewidth',3.0,'MarkerSize',12);
set(gca,'XScale','log','fontsize',30,'linewidth',3.0);
xlabel('Dose (e^-/nm^2)');
ylabel('Thickness (nm)');
set(gcf,'Position',gcf_pos);
set(gca,'Position',gca_pos);
xlim([0.8 200]);
ylim([38 51]);

%% experiment: 29 eV
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

%% Plot the deprotection images thru-z
% clc;
close all;

figure;
% errorbar(zvec',deprmean(1,simparms.zadd_npx+1:end),deprstd(1,simparms.zadd_npx+1:end),'-ob','linewidth',3.0,'MarkerSize',12);
hold on;
errorbar(zvec',deprmean(8,simparms.zadd_npx+1:end),deprstd(1,simparms.zadd_npx+1:end),'-ob','linewidth',3.0,'MarkerSize',12);
% errorbar(zvec',deprmean(6,simparms.zadd_npx+1:end),deprstd(1,simparms.zadd_npx+1:end),'-or','linewidth',3.0,'MarkerSize',12);
% errorbar(zvec',deprmean(7,simparms.zadd_npx+1:end),deprstd(1,simparms.zadd_npx+1:end),'-.b','linewidth',3.0,'MarkerSize',12);

% figure;

errorbar(zvec',deprmean(17,simparms.zadd_npx+1:end),deprstd(15,simparms.zadd_npx+1:end),'-ok','linewidth',3.0,'MarkerSize',12);
hold on;
% errorbar(zvec',deprmean(22,simparms.zadd_npx+1:end),deprstd(17,simparms.zadd_npx+1:end),'-ok','linewidth',3.0,'MarkerSize',12);
% errorbar(zvec',deprmean(18,simparms.zadd_npx+1:end),deprstd(18,simparms.zadd_npx+1:end),'-or','linewidth',3.0,'MarkerSize',12);
% errorbar(zvec',deprmean(19,simparms.zadd_npx+1:end),deprstd(19,simparms.zadd_npx+1:end),'-.b','linewidth',3.0,'MarkerSize',12);

xlim([-1 12]);
ylim([-0.1 0.8]);

set(gca,'fontsize',30,'linewidth',3.0);
box on;
