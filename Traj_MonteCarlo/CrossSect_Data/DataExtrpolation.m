%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates a new set of vibrational data by the means of
% copying an old set of data and using it as a starting point. It's good
% for incremental changes like increasing the energy range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE SURE THAT YOU"RE IN THE CORRECT DIRECTORY

%% Loading the template
clear;
load('VibrExcit_Data_Khakoo_2013.mat');
set(0,'DefaultFigureWindowStyle','Docked')

%% Make sure that the energy axes align
disp([NaN NaN Epr1.angledata(:,1)'])
disp([NaN NaN Epr2.angledata(:,1)'])
disp([NaN Epr3.angledata(:,1)'])
disp([NaN Epr4.angledata(:,1)'])
disp([NaN Epr5.angledata(:,1)'])
disp(Epr6.angledata(:,1)')


%% Turns out that the data terminates at 20 eV. Make a 3D matrix to see how things go
matDCS(:,:,6) = Epr6.angledata(:,2:end);
matDCS(:,:,:) = NaN;

matDCS(3:end,1:size(Epr1.angledata,2)-1,1) = Epr1.angledata(:,2:end);
matDCS(3:end,1:size(Epr2.angledata,2)-1,2) = Epr2.angledata(:,2:end);
matDCS(2:end,:,3) = Epr3.angledata(:,2:end);
matDCS(2:end,:,4) = Epr4.angledata(:,2:end);
matDCS(2:end,:,5) = Epr5.angledata(:,2:end);
matDCS(:,:,6) = Epr6.angledata(:,2:end);

%% See if there are trends. pick a mode
% matDCS(angle, mode, incidentEergy)
figure(7201);
for ii = 1:4
    figure(7200);
    subplot(2,2,ii);
    MODE = ii;
    [X,Y] = meshgrid([2 3 5 10 15 20],[10 15 20 30 40 50 60 70 80 90 110 130]);
    s = surf(X,Y,squeeze(matDCS(:,MODE,:)*100));
    ylabel('\theta(degrees)');xlabel('E_{incidence}(eV)');
    title(strcat('DCS of mode ',num2str(MODE)));
    s.EdgeColor = 'none';
    axis([2 20 10 130 -inf inf])
    %daspect([1 5 0]);
end

%% Try extrapolation with a power times an exponential function
% matDCS(angle, mode, incidentEergy)
CURR_ANGLE_IDX = 4;
CURR_MODE_IDX = 4;
 fitOptions =   fitoptions('Method','NonlinearLeastSquares',...
     'Lower',[0 0 0 ],...
     'Upper',[1 10 2 ],...
     'StartPoint',[max(squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:))) 0.5 0.05 ]);
 fitType    =   fittype('a*(x^b)*exp(-c*x)/((b/c)^b*exp(-b))','options',fitOptions);
 xx = [2 3 5 10 15 20]';
 yy = squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:));
 [testResult, testGOF] = fit(xx(~isnan(yy)),yy(~isnan(yy)),fitType);
 figure(7202);
 hold off
 plot(testResult);
 hold on
 plot([2 3 5 10 15 20]',squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:)))
 disp(testResult)
 
 %% Loop it to make sure that it works
 rSquaredResults = ones([12 4]);
 rSquaredResults(:,:) = NaN;
 bResults = rSquaredResults;
 cResults = bResults;
 bConfint(:,:,2) = bResults;
 cConfint(:,:,2) = bResults;
 model = 'a*(x^b)*exp(-c*x)/((b/c)^b*exp(-b))';

 xx = [2 3 5 10 15 20]';
 for CURR_ANGLE_IDX = 2:12
    for CURR_MODE_IDX = 1:4
         fitOptions =   fitoptions('Method','NonlinearLeastSquares',...
             'Lower',[0 0 0 ],...
             'Upper',[1 10 1 ],...
             'StartPoint',[max(squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:))) 0.5 0.05 ]);
          fitType    =   fittype(model,'options',fitOptions);
         yy = squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:));
        [testResult, testGOF] = fit(xx(~isnan(yy)),yy(~isnan(yy)),fitType);
        
        % if the goodness of fit is too low (mode 1 at 15 degrees) use only
        % an exponential funtion. 
        % This is actually a bad criterion. Reducing the model alwasy
        % results in worse fit
        if (testGOF.rsquare<0)
            fitOptions =   fitoptions('Method','NonlinearLeastSquares',...
             'Lower',[0 0 0 ],...
             'Upper',[1 0 1 ],...
             'StartPoint',[max(squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:))) 0 0.05 ]);
          fitType    =   fittype(model,'options',fitOptions);
            yy = squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:));
            [testResult, testGOF] = fit(xx(~isnan(yy)),yy(~isnan(yy)),fitType);
            disp('Reduced model evoked')
            disp(strcat('Angle Index : ' ,num2str(CURR_ANGLE_IDX),...
                ' Mode : ',num2str(CURR_MODE_IDX)));
        end
        
        rSquaredResults(CURR_ANGLE_IDX,CURR_MODE_IDX) = testGOF.rsquare;
        currConfInt = confint(testResult);
        
        bResults(CURR_ANGLE_IDX,CURR_MODE_IDX) = testResult.b;
        bConfint(CURR_ANGLE_IDX,CURR_MODE_IDX,:) = currConfInt(:,2);
        cResults(CURR_ANGLE_IDX,CURR_MODE_IDX) = testResult.c;
        cConfint(CURR_ANGLE_IDX,CURR_MODE_IDX,:) = currConfInt(:,3);
    end
 end
 [X,Y] = meshgrid([10 15 20 30 40 50 60 70 80 90 110 130],[1 2 3 4]);
 %% Visulas
 figure(7210);
 surf(X,Y,rSquaredResults');
 xlabel('\theta(degrees)')
 ylabel('Mode')
 zlabel('Coefficient of Determination R^2 = 1-rss/SS_{tot}')
 title({'Goodness of fit for the model', model})
 axis([10 130 1 4 0.5 1])
 colormap('jet')
 colorbar;
 
 figure(7211);
 hold off
 surf(X,Y,bResults','FaceAlpha',0);
 hold on
 lc = surf(X,Y,squeeze(bConfint(:,:,1))','FaceAlpha',0.5);
 lc.EdgeColor = 'none';
 uc = surf(X,Y,squeeze(bConfint(:,:,2))','FaceAlpha',0.5);
 uc.EdgeColor = 'none';
 xlabel('\theta(degrees)')
 ylabel('Mode')
 zlabel('b')
 title({'Value of b in the model', model})
 
 figure(7212)
 hold off
 surf(X,Y,cResults','FaceAlpha',0);
 hold on
 lc = surf(X,Y,squeeze(cConfint(:,:,1))','FaceAlpha',0.5);
 lc.EdgeColor = 'none';
 uc = surf(X,Y,squeeze(cConfint(:,:,2))','FaceAlpha',0.5);
 uc.EdgeColor = 'none';
 xlabel('\theta(degrees)')
 ylabel('Mode')
 zlabel('c')
 title({'Values of c the model', model})
 figure(7214);
 plot(bResults,cResults,'x');
 xlabel('b')
 ylabel('c')
 title({'Scatter plot between b and c in model',model})
 
 %% Fitting results are correlated. A reduced model should be in place
 bArray = reshape(bResults,1,[]);
 cArray = reshape(cResults,1,[]);
 bcFit = polyfit(bArray(~isnan(bArray)),cArray(~isnan(bArray)),1);
 
fitOptions =   fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0 0],...
    'Upper',[0.2 1/5 ],...
    'StartPoint',[0.025 1/7]);
fitType    =   fittype('a + b*x','options',fitOptions);
[bcFitResult, bcFitGOF] = fit(bArray(~isnan(bArray))',cArray(~isnan(bArray))',fitType);
figure(7215)
hold off
plot(bArray,cArray,'x');
hold on
plot (bArray,bcFitResult.a + bcFitResult.b.*bArray,'-');
xlabel('b');
ylabel('c');
legend('c vs b',strcat('linear fit c =',...
    num2str(bcFitResult.a),' + ',...
    num2str(bcFitResult.b), '*b'));
title({'Model parameters b vs c in model',model})

%% Since the uncertainty of the fitting parameters are high. Reduce the model.
% This time , the peak is fixed to 6.8540
% matDCS(angle, mode, incidentEergy)
CURR_ANGLE_IDX = 4;
CURR_MODE_IDX = 4;
fitOptions =   fitoptions('Method','NonlinearLeastSquares',...
 'Lower',[0 0],...
 'Upper',[1 10],...
 'StartPoint',[max(squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:))) 1]);
fitType    =   fittype('a*(x^b)*exp(-0.1413*b*x)/((7.0771)^b*exp(-b))','options',fitOptions);
xx = [2 3 5 10 15 20]';
yy = squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:));
[testResult, testGOF] = fit(xx(~isnan(yy)),yy(~isnan(yy)),fitType);
figure(7202);
hold off
plot(testResult);
hold on
plot([2 3 5 10 15 20]',squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:)))
disp(testResult)

%% And then we have loops again

rSquaredResults = ones([12 4]);
rSquaredResults(:,:) = NaN;
bResults = rSquaredResults;
aResults = bResults;
bConfint(:,:,2) = bResults;
aConfint(:,:,2) = aResults;
model2 = 'a*(x^b)*exp(-0.1413*b*x)/((7.0771)^b*exp(-b))';
%model2 = 'a*(x^b)*exp(-0.125*b*x)/((8)^b*exp(-b))';

xx = [2 3 5 10 15 20]';
for CURR_ANGLE_IDX = 2:12
    for CURR_MODE_IDX = 1:4
        fitOptions =   fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0 0],...
            'Upper',[1 10],...
            'StartPoint',[max(squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:))) 1]);
        fitType    =   fittype(model2,'options',fitOptions);
         yy = squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:));
        [testResult, testGOF] = fit(xx(~isnan(yy)),yy(~isnan(yy)),fitType);

        rSquaredResults(CURR_ANGLE_IDX,CURR_MODE_IDX) = testGOF.rsquare;
        currConfInt = confint(testResult);

        bResults(CURR_ANGLE_IDX,CURR_MODE_IDX) = testResult.b;
        bConfint(CURR_ANGLE_IDX,CURR_MODE_IDX,:) = currConfInt(:,2);
        aResults(CURR_ANGLE_IDX,CURR_MODE_IDX) = testResult.a;
        aConfint(CURR_ANGLE_IDX,CURR_MODE_IDX,:) = currConfInt(:,1);
    end
end

 %% Visulas
 figure(7220);
 surf(X,Y,rSquaredResults');
 xlabel('\theta(degrees)')
 ylabel('Mode')
 zlabel('Coefficient of Determination R^2 = 1-rss/SS_{tot}')
 title({'Goodness of fit for the model', model2})
 axis([10 130 1 4 0.5 1])
 colormap('jet')
 colorbar;
 
 figure(7221);
 hold off
 surf(X,Y,bResults','FaceAlpha',0);
 hold on
 lc = surf(X,Y,squeeze(bConfint(:,:,1))','FaceAlpha',0.5);
 lc.EdgeColor = 'none';
 uc = surf(X,Y,squeeze(bConfint(:,:,2))','FaceAlpha',0.5);
 uc.EdgeColor = 'none';
 xlabel('\theta(degrees)')
 ylabel('Mode')
 zlabel('b')
 title({'Value of b in the model', model2})
 
 figure(7222)
 hold off
 surf(X,Y,aResults','FaceAlpha',0);
 hold on
 lc = surf(X,Y,squeeze(aConfint(:,:,1))','FaceAlpha',0.5);
 lc.EdgeColor = 'none';
 uc = surf(X,Y,squeeze(aConfint(:,:,2))','FaceAlpha',0.5);
 uc.EdgeColor = 'none';
 xlabel('\theta(degrees)')
 ylabel('Mode')
 zlabel('c')
 title({'Values of a the model', model2})
 figure(7224);
 plot(aResults,bResults,'x');
 set(gca,'XScale','Log');
 xlabel('b')
 ylabel('c')
 title({'Scatter plot between a and b in model',model2})
 
 %% The ratio between b and c is mode dependent
 % Rerun the big model cell before running this one
 figure(7300)
 hold off;
 xScale = 0:0.5:8;
 legends = cell([1 8]);
 modeResolvedResult = cell([1 4]);
 modeResolvedGOF = cell([1 4]);
 color = {'r','k','b','g'};
 
for CURR_MODE_IDX = 1:4
    fitOptions =   fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-0.05 0],...
        'Upper',[0.1 1/5 ],...
        'StartPoint',[0 1/7]);
    fitType    =   fittype('a + b*x','options',fitOptions);
    [bcFitResult, bcFitGOF] = fit(bResults(~isnan(bResults(:,CURR_MODE_IDX)),CURR_MODE_IDX),...
        cResults(~isnan(bResults(:,CURR_MODE_IDX)),CURR_MODE_IDX),fitType);
    modeResolvedResult{CURR_MODE_IDX} = bcFitResult;
    modeResolvedGOF{CURR_MODE_IDX} = bcFitGOF;
    
    plot(bResults(~isnan(bResults(:,CURR_MODE_IDX)),CURR_MODE_IDX),...
        cResults(~isnan(bResults(:,CURR_MODE_IDX)),CURR_MODE_IDX),strcat(color{CURR_MODE_IDX},'X'));
    hold on
    plot(xScale,bcFitResult.a+xScale.*bcFitResult.b,...
        strcat(color{CURR_MODE_IDX},'-'));
    legends{2*CURR_MODE_IDX-1} = strcat('Mode: ', num2str(CURR_MODE_IDX));
    legends{2*CURR_MODE_IDX} = strcat('Mode: ', num2str(CURR_MODE_IDX),...
        'c = ',num2str(bcFitResult.a), ' + ', num2str(bcFitResult.b),'x');
end
legend(legends);
xlabel('b');
ylabel('c');
%% Anndddd, the final fit

rSquaredResultsF = ones([12 4]);
rSquaredResultsF(:,:) = NaN;
bResultsF = rSquaredResults;
aResultsF = bResults;
bConfintF(:,:,2) = bResultsF;
aConfintF(:,:,2) = aResultsF;
%model2 = 'a*(x^b)*exp(-0.1413*b*x)/((7.0771)^b*exp(-b))';


xx = [2 3 5 10 15 20]';
for CURR_ANGLE_IDX = 2:12
    for CURR_MODE_IDX = 1:4
        % The mode dependent model
        cExpression = strcat('(',num2str(modeResolvedResult{CURR_MODE_IDX}.a),...
            '+',num2str(modeResolvedResult{CURR_MODE_IDX}.b),'*b)');
        modelFinal = strcat('a*(x^b)*exp(-',...
            cExpression,...
            '*x)/((b/',...
            cExpression,...
            ')^b*exp(-b))');
        % Configure the fit
        yy = squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:));
        weights = [1 1 1 4 4 10]';
        fitOptions =   fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0 0],...
            'Upper',[1 10],...
            'StartPoint',[max(squeeze(matDCS(CURR_ANGLE_IDX,CURR_MODE_IDX,:))) 1],...
            'Weights', weights(~isnan(yy)));
        fitType    =   fittype(modelFinal,'options',fitOptions);
        
        [testResult, testGOF] = fit(xx(~isnan(yy)),yy(~isnan(yy)),fitType);

        rSquaredResultsF(CURR_ANGLE_IDX,CURR_MODE_IDX) = testGOF.rsquare;
        currConfInt = confint(testResult);

        bResultsF(CURR_ANGLE_IDX,CURR_MODE_IDX) = testResult.b;
        bConfintF(CURR_ANGLE_IDX,CURR_MODE_IDX,:) = currConfInt(:,2);
        aResultsF(CURR_ANGLE_IDX,CURR_MODE_IDX) = testResult.a;
        aConfintF(CURR_ANGLE_IDX,CURR_MODE_IDX,:) = currConfInt(:,1);
    end
end
 %% Visulas
 modelF = 'a*(x^b)*exp(-c(m,b)*x)/((b/c(m,b))^b*exp(-b))';
 figure(7230);
 surf(X,Y,rSquaredResultsF');
 xlabel('\theta(degrees)')
 ylabel('Mode')
 zlabel('Coefficient of Determination R^2 = 1-rss/SS_{tot}')
 title({'Goodness of fit for the model', modelF})
 axis([10 130 1 4 0.5 1])
 colormap('jet')
 colorbar;
 
 figure(7231);
 hold off
 surf(X,Y,bResultsF','FaceAlpha',0);
 hold on
 lc = surf(X,Y,squeeze(bConfintF(:,:,1))','FaceAlpha',0.5);
 lc.EdgeColor = 'none';
 uc = surf(X,Y,squeeze(bConfintF(:,:,2))','FaceAlpha',0.5);
 uc.EdgeColor = 'none';
 xlabel('\theta(degrees)')
 ylabel('Mode')
 zlabel('b')
 title({'Value of b in the model', modelF})
 
 figure(7232)
 hold off
 surf(X,Y,aResultsF','FaceAlpha',0);
 hold on
 lc = surf(X,Y,squeeze(aConfintF(:,:,1))','FaceAlpha',0.5);
 lc.EdgeColor = 'none';
 uc = surf(X,Y,squeeze(aConfintF(:,:,2))','FaceAlpha',0.5);
 uc.EdgeColor = 'none';
 xlabel('\theta(degrees)')
 ylabel('Mode')
 zlabel('c')
 title({'Values of a the model', model2F})
 figure(7233);
 plot(aResultsF,bResultsF,'x');
 set(gca,'XScale','Log');
 xlabel('a')
 ylabel('b')
 title({'Scatter plot between a and b in model',model2F})
 
 %% The extrapolation
 % This rearrangement allows for vectorization in the cross-section
 % calculation
 C_MODE_RESOLVED = zeros([2 4]);
 for mode_i = 1:4 
     C_MODE_RESOLVED(:,mode_i) = [ modeResolvedResult{mode_i}.a;...
         modeResolvedResult{mode_i}.b];
 end
 A_MODE_ANGLE_RESOLVED = aResultsF;
 B_MODE_ANGLE_RESOLVED = bResultsF;
 
 save('DCS_models','A_MODE_ANGLE_RESOLVED','B_MODE_ANGLE_RESOLVED',...
'C_MODE_RESOLVED','modeResolvedResult')
 
 param.C_MODE_RESOLVED          =   C_MODE_RESOLVED;
 param.B_MODE_ANGLE_RESOLVED    =   B_MODE_ANGLE_RESOLVED;
 param.A_MODE_ANGLE_RESOLVED    =   A_MODE_ANGLE_RESOLVED;
 
 % Test the function
 % matDCS(angle, mode, incidentEergy)
 [angle_X,mode_Y,eInc_Z] = meshgrid(1:1:12,...
     [1 2 3 4],...
     1:1:6);
 % Vectorization doesn't work well with multidimensional arrays
 modelResults = matDCS;
 for mode_idx = 1:4
     for angle_idx = 1:12
        modelResults(angle_idx,mode_idx,:)...
            = xSecFunction(mode_idx,angle_idx,[2 3 5 10 15 20],param);
     end
 end
 
 %% Inspecting the results
 figure(7300);
 hold off;
 for idx = 1:4
     legends_verification{idx} = num2str(idx);
    loglog(reshape(matDCS(:,idx,:),[],1),...
        reshape(modelResults(:,idx,:),[],1),'X');
    hold on;
 end
 title({'Model Deviation from Measured Values','Segment by mode'});
 xlabel('Measured DCS');
 ylabel('Model DCS');
 legend(legends_verification);
 
 %% Extend the numbers
extendedResults = zeros([12 4 3]);
for mode_idx = 1:4
     for angle_idx = 1:12
        extendedResults(angle_idx,mode_idx,:)...
            = xSecFunction(mode_idx,angle_idx,[30 40 50],param);
     end
end

%% Putting it next to the data and see how well it does
% matDCS(angle, mode, incidentEergy)
figure(7400);
for ii = 1:4
    figure(7400);
    subplot(2,2,ii);
    MODE = ii;
    
    hold off
    % The original data
    [X,Y] = meshgrid([2 3 5 10 15 20],[10 15 20 30 40 50 60 70 80 90 110 130]);
    s1 = surf(X,Y,squeeze(matDCS(:,MODE,:)*100));
    hold on
    s1.EdgeColor = 'none';
    
    % The extended data
    [X,Y] = meshgrid([30 40 50],[10 15 20 30 40 50 60 70 80 90 110 130]);
    s2 = surf(X,Y,squeeze(extendedResults(:,MODE,:)*100));
    s2.EdgeColor = 'none';
    
    
    ylabel('\theta(degrees)');xlabel('E_{incidence}(eV)');
    title(strcat('DCS of mode ',num2str(MODE)));
    set(gca,'ZScale','Log');
    
    axis([2 50 10 130 -inf inf])
    %daspect([1 5 0]);
end

%% The ICS. I will just extrapolate from the ics data.
% The same model should work
figure(7500);
[X,Y] = meshgrid([2 3 5 10 15 20],[1 2 3 4]);
surf(X,Y,ics(:,2:end)');
set(gca,'ZScale','Log');

%% Just use the reduced model and get something
xx = [2 3 5 10 15 20]';

rSquaredResultsICS = ones([1 4]);
rSquaredResultsICS(:) = NaN;
bResultsICS = rSquaredResultsICS;
aResultsICS = bResultsICS;
bConfintICS(:,2) = bResultsICS;
aConfintICS(:,2) = aResultsICS;

for CURR_MODE_IDX = 1:4
    % The mode dependent model
    cExpression = strcat('(',num2str(modeResolvedResult{CURR_MODE_IDX}.a),...
        '+',num2str(modeResolvedResult{CURR_MODE_IDX}.b),'*b)');
    modelFinal = strcat('a*(x^b)*exp(-',...
        cExpression,...
        '*x)/((b/',...
        cExpression,...
        ')^b*exp(-b))');
    % Configure the fit
    yy = ics(1:6,CURR_MODE_IDX+1);
    weights = [0.25 0.25 0.25 4 4 100]';
    fitOptions =   fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[0 0],...
        'Upper',[5 10],...
        'StartPoint',[max(ics(:,CURR_MODE_IDX+1)) 1],...
        'Weights',weights(~isnan(yy)));
    fitType    =   fittype(modelFinal,'options',fitOptions);
    [testResult, testGOF] = fit(xx(~isnan(yy)),yy(~isnan(yy)),fitType);

    rSquaredResultsICS(CURR_MODE_IDX) = testGOF.rsquare;
    currConfIntICS = confint(testResult);

    bResultsICS(CURR_MODE_IDX) = testResult.b;
    bConfintICS(CURR_MODE_IDX,:) = currConfInt(:,2);
    aResultsICS(CURR_MODE_IDX) = testResult.a;
    aConfintICS(CURR_MODE_IDX,:) = currConfInt(:,1);
    
    figure(7501);
    subplot(2,2,CURR_MODE_IDX);
    hold off
    plot(xx,yy);
    hold on
    param.C_MODE_RESOLVED          =   C_MODE_RESOLVED;
    param.B_MODE_ANGLE_RESOLVED    =   bResultsICS;
    param.A_MODE_ANGLE_RESOLVED    =   aResultsICS;
    plot(0:1:50,xSecFunction(CURR_MODE_IDX,1,0:1:50,param));
    title(strcat('ICS of mode ',num2str(CURR_MODE_IDX)));
    xlabel('Energy (eV)');
    ylabel('ICS A^{-2}');
    
end

%% Update the ICS array
ics(7:9,1) = [30 40 50];
for ii = 2:5
    ics(7:9,ii) = xSecFunction(ii-1,1,[30 40 50],param);
end
save('Khakoo_2013_50eV.mat','ics');

%% Generate new data entries
m = matfile('Vibr_Khakoo_2013_x50','Writable',true);
m.ics = ics;
angleArray = [10 15 20 30 40 50 60 70 80 90 110 130];


%% 30 eV
Epr7=Epr6;
Epr7.Eo = 30;
Epr7.angledata = zeros(11,5);

Epr7.angledata(:,1) = angleArray(2:end);
Epr7.angledata(:,2:end) = squeeze(extendedResults(2:end,:,1));

%% 40 eV
Epr8=Epr6;
Epr8.Eo = 40;
Epr8.angledata = zeros(11,5);

Epr8.angledata(:,1) = angleArray(2:end);
Epr8.angledata(:,2:end) = squeeze(extendedResults(2:end,:,2));

%% 50 eV
Epr9=Epr6;
Epr9.Eo = 50;
Epr9.angledata = zeros(11,5);

Epr9.angledata(:,1) = angleArray(2:end);
Epr9.angledata(:,2:end) = squeeze(extendedResults(2:end,:,3));

%% Save the results
m.Epr7 = Epr7;
m.Epr8 = Epr8;
m.Epr9 = Epr9;

