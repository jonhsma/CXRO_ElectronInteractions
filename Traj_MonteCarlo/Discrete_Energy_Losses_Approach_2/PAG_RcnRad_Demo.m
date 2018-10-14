%% PAG reaction radius demo
clc;
close all;

x=-5:0.5:5;
[x,y,z]=meshgrid(x);

rcnrad=2;
rcnsphere=zeros(size(z));

z(z>=-sqrt(rcnrad^2-(x.^2+y.^2)) & z<=sqrt(rcnrad^2-(x.^2+y.^2)))=1;

tmpidx=find(z==1);
[xidx,yidx,zidx]=ind2sub(size(rcnsphere),tmpidx);

%%
% r=1;
% phi=linspace(0,pi,30);
% theta=linspace(0,2*pi,40);
% [phi,theta]=meshgrid(phi,theta);
% 
% x=r*sin(phi).*cos(theta);
% y=r*sin(phi).*sin(theta);
% z=r*cos(phi); 