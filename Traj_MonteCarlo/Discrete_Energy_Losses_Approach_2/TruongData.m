% PMMA:
E_pmma = [20   30  40  40  50  70  70 90  120 200 300 400 500 700 1000 10000] % eV
Sp = [0.1 0.3 0.6 0.8 1.4 1.8 2  2.6 3.0 3.2 3.0 2.4 2.2 1.6 1.38 0.2]; % eV/Anstrom
Sp_pmma=Sp.*10; % eV/nm
Sp_pmma_unit='eV/nm';

% PolyStyrene
E_ps = [20   30  45  60  75  100 160 200 280 400 700 1000 2000 10000]; % eV
Sp = [0.1 0.3 0.6 1.0 2.0 2.8 3.0 3.0 2.8 2.2 1.4 1.1  0.8  0.1]; % eV/Angstrom
Sp_ps=Sp.*10;
Sp_ps_unit='eV/nm';

figure;
plot(E_pmma,Sp_pmma,'-ob');
hold on; plot(E_ps,Sp_ps,'-or');
set(gca,'XScale','log');