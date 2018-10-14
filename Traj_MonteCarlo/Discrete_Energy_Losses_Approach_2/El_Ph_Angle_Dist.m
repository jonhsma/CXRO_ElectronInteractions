function El_Ph_Angle_Dist(scattdata,ntrials,E)

eps0=scattdata.eps0;
epsinf=scattdata.epsinf;
E_optphonon=scattdata.hbarw;

theta=[];
for i = 1:ntrials
    [imfp_T,imfp_creation,imfp_annihilation,theta(i)]=scattdata.imfp_func(E,eps0,epsinf,E_optphonon);
end

[n,x]=hist(theta,100);

figure;
bar(x,n./trapz(x,n),1);
xlabel('Theta (degrees)');
ylabel('Normalized Count');