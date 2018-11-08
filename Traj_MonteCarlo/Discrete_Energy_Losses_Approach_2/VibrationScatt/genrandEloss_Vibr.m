function outdata=genrandEloss_Vibr(scattdata,incidentE)

E   =   scattdata.ics(:,1);
ics =   scattdata.ics(:,2:end);

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
%icsvec=[];
if ~isempty(idx2) % found exact match for Eo
    icsvec=ics(idx2,:);
    Eochoose=idx2;
else % did not find an exact match for Eo
    vec1=ics(idx1,:);
    vec2=ics(idx3,:);
    icsvec = zeros([1 size(vec1,2)]);
    for count=1:size(vec1,2)
        p=polyfit([E(idx1) E(idx3)],[vec1(1,count) vec2(1,count)],1);
        icsvec(1,count)=polyval(p,incidentE);
    end
    Eochoose=[idx1 idx3];
end

%% Monte-Carlo the mode of vibration involved
%%% I'm not convinced with this implementation of monte carlo
%{
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
%}

%%% Generate the mode culmulative probability distribution
activeModes     =   find(~isnan(icsvec));   % Holds the active modes
activeICS       =   icsvec(activeModes);    % Holds the cross-section active modes
modeCumICS_CDF  =   cumsum(activeICS);      % The cumulativeICS of the active modes only
modeCumICS_CDF  =   modeCumICS_CDF/modeCumICS_CDF(end);
%%%
modeIdx     =   activeModes(find(modeCumICS_CDF>rand,1));

Elchoose    =   modeIdx; % the EL value [which vibrational mode?]
Eloss       =   (scattdata.El_min(Elchoose)+...
    scattdata.El_max(Elchoose))/2;

%%%% now pick the scattering angle
dcsvec=[];
%%% theta_int is theta scale in dgree
theta_int=10:130;
if length(Eochoose)==1 % found exact match for Eo
    thetavec=eval(sprintf('scattdata.Epr%d.angledata(:,1)',Eochoose));
    dcsvec=eval(sprintf('scattdata.Epr%d.angledata(:,Elchoose+1)',Eochoose));
    dcsvec=interp1(thetavec,dcsvec,theta_int,'spline','extrap');
else % did not find an exact match for Eo
    theta1=eval(sprintf('scattdata.Epr%d.angledata(:,1)',Eochoose(1)));
    theta2=eval(sprintf('scattdata.Epr%d.angledata(:,1)',Eochoose(2)));
    
    vec1=eval(sprintf('scattdata.Epr%d.angledata(:,Elchoose+1)',Eochoose(1))); 
    vec1=interp1(theta1,vec1,theta_int,'spline','extrap');
%     vec1=vec1'; % turn it into a row vector for the loop below to work
    
    vec2=eval(sprintf('scattdata.Epr%d.angledata(:,Elchoose+1)',Eochoose(2)));
    vec2=interp1(theta2,vec2,theta_int,'spline','extrap');
%     vec2=vec2';
    
    for count=1:size(vec1,2)
        p=polyfit([E(idx1) E(idx3)],[vec1(1,count) vec2(1,count)],1);
        dcsvec(1,count)=polyval(p,incidentE);
    end
end

dcs_cdf     =   cumsum(sin(theta_int/180*pi).*dcsvec);
dcs_cdf     =   dcs_cdf-dcs_cdf(1);
dcs_cdf     =   dcs_cdf./dcs_cdf(end);

theta_rand=randgen(theta_int,dcs_cdf,1);
%%% theta_rand is in degree here. I have not idea what is happening.
%theta_rand=pi-theta_rand; % data is deflection, so pi-() gives you elevation angle
phi_rand=2*pi*rand;

%%%%% finally, calculate the total ics;
Etmp=scattdata.ics(:,1);
ics2=scattdata.ics(:,2:end);
ics2(isnan(ics2))=0;
ics2=sum(ics2,2);
if length(Eochoose)==1
    icsval=ics2(Eochoose);
else
    p=polyfit([Etmp(Eochoose(1)) Etmp(Eochoose(2))],[ics2(Eochoose(1)) ics2(Eochoose(2))],1);
    icsval=polyval(p,incidentE);
end

outdata.Eloss=Eloss;
outdata.ics=icsval*scattdata.ics_mult;
outdata.theta=theta_rand.*pi/180;
outdata.phi=phi_rand;

% flag=1;