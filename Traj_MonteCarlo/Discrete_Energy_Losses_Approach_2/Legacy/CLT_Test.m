%% CLT test
clear;
clc;
close all;

n=100000;

for i = 1:n
    randvec(i)=mean(poissrnd(2,1,100));
end

std(randvec)
std(poissrnd(2,1,1000))/sqrt(n)

figure;
hist(randvec,50);

figure
hist(poissrnd(2,1,1000),100);