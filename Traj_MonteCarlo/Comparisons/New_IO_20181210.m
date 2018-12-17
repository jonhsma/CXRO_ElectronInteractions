%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This comparison does a few things
%   1   Test the new loading capabilities of the comparison script
%   2   Make sure that things are consistent from Dec 5 to Dec 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Making sure that things load
loadTest.runsToCompare =...
    {'20181206_Ref'};
loadTest.incEnergy =...
    {'80'};
loadTest.doseStr  = {'0.00'};
loadTest.nTrials  = {1100,1100};
loadTest.res      =   0.5;
loadTest.cuOption =   'cumcount';
loadTest.distOption =   'pdf';
loadTest.rangeLimit =   10;

loadTestResult = Comparison_ScatteringSim_Acid(loadTest);

%% Analyzing the run from Dec 6
loadDec6.runsToCompare =...
    {'20181206_Ref';'20181206_Ref';'20181206_Ref';...
    '20181206_Ref';'20181206_Ref';'20181206_Ref'};
loadDec6.incEnergy =...
    {'20','30','40','50','65','80'};
loadDec6.doseStr  = {'0.00';'0.00';'0.00';'0.00';'0.00';'0.00'};
loadDec6.nTrials  = {1100,1100,1100,1100,1100,1100};
loadDec6.res      =   0.5;
loadDec6.cuOption =   'cumcount';
loadDec6.distOption =   'pdf';
loadDec6.rangeLimit =   10;

dec6RefResult = Comparison_ScatteringSim_Acid(loadDec6);

%% Comparing the runs on Dec5 and Dec6
% Get the 1205 data
ContingentScriptFor20181205
tabulateSummary(resultObject)
tabulateSummary(dec6RefResult)

%% Some closer look at the Dec 06 data
figure(9004);
hold off
histogram(dec6RefResult{1}.rabs,30,'Normalization','pdf')
hold on
histogram(dec6RefResult{2}.rabs,30,'Normalization','pdf')
legend('20eV','30eV')
fprintf('\n%15s%10s%10s','Incident Energy','RMS(r)','<r>');
fprintf('%15s%10.2d%10.2d\n','20eV',...
    mean(dec6RefResult{1}.rabs.^2).^0.5,...
    mean(dec6RefResult{1}.rabs))
fprintf('%15s%10.2d%10.2d\n','30eV',...
    mean(dec6RefResult{2}.rabs.^2).^0.5,...
    mean(dec6RefResult{2}.rabs))