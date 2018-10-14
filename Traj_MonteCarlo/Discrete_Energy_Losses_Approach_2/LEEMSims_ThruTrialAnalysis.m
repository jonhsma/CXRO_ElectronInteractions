%% LEEM Simdata thru-trial analysis from data saved in .dat files: Load the data
% clear;
clc;
close all;

datamat=[   29		2		48.2542    	1.0934   	50.0000     0
            29		4		45.9005    	0.4647   	49.8862    	0.3599
            29		8		44.5565    	0.2677   	48.3596    	1.5223
            29		10		43.8046   	0			46.1515		0
            29		12		43.5577    	0.4421   	46.1063    	0.5267
            29		16		42.8225    	0.4137   	45.2511    	0.4814
            29		32		41.5717   	0			43.6337		0
            29		64		40.2402   	0			42.3431		0
            49		2		47.3347   	0			50.0000		0
            49		4		44.9703   	0			47.6733		0
            49		8		44.1962   	0			46.8090		0
            49		16		41.5624   	0			43.8564		0
            49		32		40.7913   	0			42.9131		0
            49		64		39.9086   	0			42.1600		0
            91		2		44.7288    	0.2925   	48.7583    	1.3124
            91		4		43.2263    	0.5442   	45.6158    	0.5826
            91		8		41.8334    	0.2282   	43.8991    	0.2411
            91		10		41.5523   	0			43.8713		0
            91		16		40.6565    	0.5238   	42.7582    	0.3995
            91		32		39.6780   	0			41.6381		0
            91		64		37.5925   	0			39.5539		0];
        
%% create a structure "rtdata" from the "datamat" matrix
clc;
close all;

rtdata={};

E=datamat(:,1);

Eval=29;
idx=find(E==Eval);
datalen=length(rtdata);
rtdata{datalen+1}.dose=datamat(idx,2);
rtdata{datalen+1}.rtmean_0p3=datamat(idx,3);
rtdata{datalen+1}.rtstd_0p3=datamat(idx,4);
rtdata{datalen+1}.rtmean_0p5=datamat(idx,5);
rtdata{datalen+1}.rtstd_0p5=datamat(idx,6);
rtdata{datalen+1}.legstr=sprintf('E=%.1f',Eval);

Eval=49;
idx=find(E==Eval);
datalen=length(rtdata);
rtdata{datalen+1}.dose=datamat(idx,2);
rtdata{datalen+1}.rtmean_0p3=datamat(idx,3);
rtdata{datalen+1}.rtstd_0p3=datamat(idx,4);
rtdata{datalen+1}.rtmean_0p5=datamat(idx,5);
rtdata{datalen+1}.rtstd_0p5=datamat(idx,6);
rtdata{datalen+1}.legstr=sprintf('E=%.1f',Eval);

Eval=91;
idx=find(E==Eval);
datalen=length(rtdata);
rtdata{datalen+1}.dose=datamat(idx,2);
rtdata{datalen+1}.rtmean_0p3=datamat(idx,3);
rtdata{datalen+1}.rtstd_0p3=datamat(idx,4);
rtdata{datalen+1}.rtmean_0p5=datamat(idx,5);
rtdata{datalen+1}.rtstd_0p5=datamat(idx,6);
rtdata{datalen+1}.legstr=sprintf('E=%.1f',Eval);

%% Plot values in "rtdata" structure
clc;
close all;

pltstyles={'-ob','-ok','-or','-og','--ob','--ok','--or','--og'};
p_0p3=[];
p_0p5=[];

gca_pos=[0.1300    0.1100    0.7750    0.8150];
gcf_pos=[ 680         181        1076         797];

%%%% threshold=0.3
figure;hold on;xlabel('Dose (e^-/nm^2)');ylabel('Thickness (nm)');
set(gca,'fontsize',30,'linewidth',3.0,'XScale','log');
set(gcf,'Position',gcf_pos);
set(gca,'Position',gca_pos);

for i = 1:length(rtdata)
    dosevec=rtdata{i}.dose;
    rtvec=rtdata{i}.rtmean_0p3;
    rt_lb=0;    rt_ub=49;
    p=polyfit(log10(dosevec(rtvec<=rt_ub&rtvec>=rt_lb)),rtvec(rtvec<=rt_ub&rtvec>=rt_lb),1);
    p_0p3(i)=p(1);
    
    errorbar(rtdata{i}.dose,rtdata{i}.rtmean_0p3,rtdata{i}.rtstd_0p3,pltstyles{i},'linewidth',3.0,'markersize',12);
    drawnow;
end

%%%% threshold=0.5
figure;hold on;xlabel('Dose (e^-/nm^2)');ylabel('Thickness (nm)');
set(gca,'fontsize',30,'linewidth',3.0,'XScale','log');
set(gcf,'Position',gcf_pos);
set(gca,'Position',gca_pos);

for i = 1:length(rtdata)
    dosevec=rtdata{i}.dose;
    rtvec=rtdata{i}.rtmean_0p5;
    rt_lb=43;    rt_ub=49;
    dosefit=dosevec(rtvec<=rt_ub&rtvec>=rt_lb);
    rtfit=rtvec(rtvec<=rt_ub&rtvec>=rt_lb);
    p=polyfit(log10(dosefit),rtfit,1);
    p_0p5(i)=p(1);
    
    errorbar(rtdata{i}.dose,rtdata{i}.rtmean_0p5,rtdata{i}.rtstd_0p5,pltstyles{i},'linewidth',3.0,'markersize',12);
    drawnow;
end

[p_0p3;p_0p5]