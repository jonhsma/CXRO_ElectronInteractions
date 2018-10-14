%% testbench to figure out x,y,z and angles
clear;
clc;
close all;


x=0;
y=0;
z=0;

angle(1)=3*pi/4;

nsteps=100;
pathlen=10;
phi=pi/4;

for i = 2:nsteps
    phi=2*pi*rand;
%     phi=0;
    angle(i)=angle(i-1)+pi/2*rand;
    z(i)=z(i-1)+pathlen*cos(angle(i));
    x(i)=x(i-1)+pathlen*sin(angle(i))*cos(phi);
    y(i)=y(i-1)+pathlen*sin(angle(i))*sin(phi);
end

figure;
plot3(x,y,z,'o-b')
xlabel('x');
ylabel('y');
zlabel('z');
grid on;