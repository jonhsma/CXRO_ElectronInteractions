%%%% Some un organized scratch work
%% Dec 28
[escapeEvents,incidentEvents] = extractEscapeEvents(energyScanArchive);
coords = zeros([6575,3]);
%{
for ii = 1:6575
    this = escapeEvents(ii);coords(ii,:)=[cos(this.theta_in),sin(this.theta_in)*cos(this.phi_in),sin(this.theta_in)*sin(this.phi_in)];
end
%}
for ii = 1:6575
    this = escapeEvents(ii);coords(ii,:)=[cos(this.theta_out),sin(this.theta_out)*cos(this.phi_out),sin(this.theta_out)*sin(this.phi_out)];
end

figure(3002);plot3(coords(:,3),coords(:,2),coords(:,1),'.')
daspect([1 1 1])
axis([-1 1 -1 1 -1 1])
title('Direction of escaped electrons projected onto a unit sphere')

energyAtEscape = zeros([1 6575]);
for ii = 1:6575
this = escapeEvents(ii);energyAtEscape(ii) = this.Ein;
end
figure(3003); histogram(energyAtEscape)
title({'Simulate Photoemission spectrum';'Primary KE = 85 eV';'primaries from the first 5 nm considered'})
xlabel('KE(eV)')
ylabel('Counts')
title({'Simulated Photoemission spectrum';'Primary KE = 85 eV';'primaries from the first 5 nm considered'})

figure(3013); histogram(energyAtEscape)
title({'Simulate Photoemission spectrum';'Primary KE = 85 eV';'primaries from the first 5 nm considered'})
xlabel('KE(eV)')
ylabel('Counts')
title({'Simulated Photoemission spectrum';'Primary KE = 85 eV';'primaries from the first 5 nm considered'})
set(gca, 'YScale', 'log')
%% Dec 29
figure(3004);
e0 = 80;
hold off
plot3(coords(energyAtEscape<e0,3),coords(energyAtEscape<e0,2),coords(energyAtEscape<e0,1),'.')
hold on
plot3(coords(energyAtEscape>=e0,3),coords(energyAtEscape>=e0,2),coords(energyAtEscape>=e0,1),'.')
daspect([1 1 1])
axis([-1 1 -1 1 0 1])
title('Direction of escaped electrons projected onto a unit sphere')

escapeTheta = zeros([1 6575]);
figure(3005);
hold off
hLow = histogram([escapeEvents(energyAtEscape<e0).theta_out],'Normalization','pdf');
hold on 
hHigh = histogram([escapeEvents(energyAtEscape>=e0).theta_out],'Normalization','pdf');
plot(0:pi/400:pi/2,sin((0:pi/400:pi/2)));
title('Theta distribution of escaped electrons')
legend(strcat('< ',num2str(e0),' eV emissions'),...
    strcat('>= ',num2str(e0),' eV emissions'),...
    'Random (sin(\theta)) reference');
axis([0 pi/2 -inf inf])
xlabel('\theta');
ylabel('Relative probability per unit \theta');

figure(3006);
hold off
barcenterLow = (hLow.BinEdges(1:end-1)+hLow.BinEdges(2:end))/2;
plot(barcenterLow,hLow.Values./barcenterLow)
hold on
barcenterHigh = (hHigh.BinEdges(1:end-1)+hHigh.BinEdges(2:end))/2;
plot(barcenterHigh,hHigh.Values./barcenterHigh)
legend(strcat('< ',num2str(e0),' eV emissions'),...
    strcat('>= ',num2str(e0),' eV emissions'));
axis([0 pi/2 -inf inf])
title('Probability distribution per unit solid angle')
xlabel('\theta');
ylabel('Relative probability per unit solid angle');

%% Another sanity check (Dec 29) what is the depth distribution of the escaped electrons
e1=80;

zIdx = 1:1:6575;
zIdx = zIdx*3;
xyzArray = [incidentEvents.xyz_init];
zArray = xyzArray(zIdx);

figure(4001);
hold off
histogram(zArray(energyAtEscape<e1));
hold on
histogram(zArray(energyAtEscape>=e1));
title({'Distributionof depth of photoemission origin';...
    'for all escaped electrons'});
xlabel('z(nm) Resist exists below 0');
ylabel('raw counts')
legend(strcat('KE <',num2str(e1),'eV'),...
    strcat('KE >=',num2str(e1),'eV'),...
    'Location','northwest');


figure(4002);
hold off
hDLow = histogram(zArray(energyAtEscape<e1));
hold on
histogram(zArray(energyAtEscape>=e1));
title({'Distributionof depth of photoemission origin';...
    'for all escaped electrons'});
xlabel('z(nm) Resist exists below 0');
ylabel('raw counts')
set(gca, 'YScale', 'log')
legend(strcat('KE <',num2str(e1),'eV'),...
    strcat('KE >=',num2str(e1),'eV'),...
    'Location','northwest');

%% Get the log slope of the depth distribution of escaped secondaries
barcenterDLow = (hDLow.BinEdges(1:end-1)+hDLow.BinEdges(2:end))/2;
f = fit (barcenterDLow(barcenterDLow<-1)',...
    hDLow.Values(barcenterDLow<-1)','exp1');
figure(5200);
hold off
bar(barcenterDLow,hDLow.Values);
hold on;
plot(barcenterDLow(barcenterDLow<-1),...
    f.a.*exp(f.b.*barcenterDLow(barcenterDLow<-1)));
title({'Depth distribution of origin';strcat('KE <',num2str(e1),' eV');...
    strcat('\lambda = ',num2str(1/f.b),' nm')})
xlabel('z(nm) Resist exists below 0');
ylabel('raw counts')
legend({'Depth distribution';...
    strcat('Fit: e^{z/\lambda};\lambda = ',num2str(1/f.b))});

figure(5201)
hold off
hE = histogram(energyAtEscape);
barCenterE = (hE.BinEdges(1:end-1)+hE.BinEdges(2:end))/2;
f = fit (barCenterE(barCenterE<=40&barCenterE>=17)',...
    hE.Values(barCenterE<=40&barCenterE>=17)','exp1');
hold on
plot(barCenterE(barCenterE<=40&barCenterE>=17),...
    f.a.*exp(barCenterE(barCenterE<=40&barCenterE>=17)*f.b));
title({'Energy distribution of escaped electrons';...
    strcat('\lambda_E = ',num2str(abs(1./f.b)),'eV')})
xlabel('KE (eV)')
legend('Simulated energy spectrum','Local exponential fit')



