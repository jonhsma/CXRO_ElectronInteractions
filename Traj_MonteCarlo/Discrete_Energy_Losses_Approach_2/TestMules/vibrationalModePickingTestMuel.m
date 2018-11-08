%%% This script/function is a largely abridged verion of the vibrational scattering engine
%%% It verifies that the correct mode is selected using integrated
%%% scattering cross-setion as the probability weight
%%% One can specify the variable name or use the variable vibrScattdata in
%%% the workspace

function outdata=vibrationalModePickingTestMuel(varargin)

addpath('..\');
N_TRIALS        =   720;
FIG_NUMBER   =   3200;
%%% Check if scattering data is there
switch(nargin)
    case 0
        if ~exist(vibrScattdata,'var')
            disp('No vibrational scattering data in workspace or input');
            return
        else
            incidentE = rand*50+30;
        end            
    case 1
        if ~exist(vibrScattdata,'var')
            disp('No vibrational scattering data in workspace or input');
        else
            incidentE = varargin{1};
        end
    case 2
        incidentE       =   varargin{1};
        vibrScattdata   =   varargin{2};
    case 3
        incidentE       =   varargin{1};
        vibrScattdata   =   varargin{2};
        N_TRIALS        =   varargin{3};
    otherwise
        incidentE       =   varargin{1};
        vibrScattdata   =   varargin{2};
        N_TRIALS        =   varargin{3};
        FIG_NUMBER      =   varargin{4};
end
E   =   vibrScattdata.ics(:,1);
ics =   vibrScattdata.ics(:,2:end);

%%% check if there is ics

%{
% Efit=linspace(0,max(E),1000);
% icsfit=zeros([length(Efit) size(ics,2)]);
% 
% for i = 1:size(ics,2)
%     icstmp=ics(:,i);
%     icsfit(:,i)=interp1(E(~isnan(icstmp)),icstmp(~isnan(icstmp)),Efit,'linear','extrap')';
% end
% 
% E=Efit;
% ics=icsfit;
%}
%% Find the closest energy entry in the crosssection table
idx1    =   find(E<incidentE);
[eJustBelow,idxtmp]=min(abs(E(idx1)-incidentE));
idx1=idx1(idxtmp);

idx2=find(E==incidentE);

idx3=find(E>incidentE);
[eJustAbove,idxtmp]=min(abs(E(idx3)-incidentE));
idx3=idx3(idxtmp);

%% generate the ICS vector at desired energy
icsvec=[];
if ~isempty(idx2) % found exact match for Eo
    icsvec=ics(idx2,:);
    Eochoose=idx2;
else % did not find an exact match for Eo
    vec1=ics(idx1,:);
    vec2=ics(idx3,:);
    icsvec = zeros([1 size(vec1,2)]);
    for count=1:size(vec1,2)
        p=polyfit([eJustBelow eJustAbove],[vec1(1,count) vec2(1,count)],1);
        icsvec(1,count)=polyval(p,incidentE);
    end
    Eochoose=[idx1 idx3];
end

%% Monte-Carlo test for the mode of vibration involved

resOriginal     =   zeros([1 N_TRIALS]);
resNew          =   zeros([1 N_TRIALS]);

for ii = 1: N_TRIALS
    ics2        =   icsvec./sum(icsvec(~isnan(icsvec)));
    randval     =   rand; % pick a uniformly distributed random number between 0 and 1
    idx         =   find(ics2<=randval);
    ics3        =   ics2;
    ics3(setdiff([1:length(ics3)],idx))=NaN;
    [minval,idxtmp]=min(abs(ics3-randval));
    if length(idxtmp)>1
        randval2=randi([1 length(idxtmp)]);
        idxtmp=idxtmp(randval2);
    end
    
    resOriginal(ii)     =   idxtmp;


    %%% Generate the mode culmulative probability distribution
    modeCumICS_CDF  =   cumsum(icsvec);
    modeCumICS_CDF  =   modeCumICS_CDF/modeCumICS_CDF(end);

    resNew(ii)      =   find(modeCumICS_CDF>rand,1);    
end
%% Output
figure(FIG_NUMBER);
hold off;
plot(icsvec/sum(icsvec),'r-X');
hold on
title({'Probability distributions of vibrational modes',...
    strcat('Incident energy=',num2str(incidentE),'eV')})
histogram(resOriginal,'Normalization','pdf');
histogram(resNew,'Normalization','pdf');
legend('ICS Data','Original','New')

