function resultObject = LoadArchives(runsToCompare,nRuns,incEngy,nTrials,doseStr)
  %% Pull up the acid distribution
    resultObject = cell([1 nRuns]);
    currentWorkSpace = '';
    for runCnt = 1:size(runsToCompare,1)
        counter = 1;
        counter_f = 1;
        resultObject{runCnt}.acid_xyz_accul=[];
        resultObject{runCnt}.acid_fine_xyz_accul=[];

        fprintf('Currently Processing')
        disp(strcat(runsToCompare(runCnt,:),...
            '; Energy :',incEngy{runCnt}));
        
        %% The file names in different formats
        % The parallelsim enabled ones
        fileName_Para   =   strcat('..\..\..\..\',...
                'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\',...
                runsToCompare{runCnt},'\ScanArchive.mat');
        % The energy separated ones
        fileName_Sep    =   strcat('..\..\..\..\',...
                'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\',...
                runsToCompare{runCnt},...
                '\ScanArchive_',incEngy{runCnt},'.mat');
        
        if exist(fileName_Para,'file')
            %% Scan archives for parallelism enabled runs
            if ~strcmp(currentWorkSpace,runsToCompare{runCnt})
                load(fileName_Para,'\ScanArchive.mat');
                currentWorkSpace = runsToCompare{runCnt};
            end
            if ~exist('scanArchive','var')
                disp('Scan Archive exists but not accessible');
                continue
            end
            %%% Find the right incident energy
            tic
            for ii = 1:size(scanArchive,2)
                if scanArchive{ii}{1}.energy == str2double(incEngy{runCnt})
                    for jj = 1:min(size(scanArchive{ii},2),nTrials{runCnt})
                        acid_xyz        =   scanArchive{ii}{jj}.acid_xyz;
                        acid_fine_xyz   =   scanArchive{ii}{jj}.activation_xyz;
                        if size(acid_xyz,1)>=1
                            resultObject{runCnt}.acid_xyz_accul(counter:counter + size(acid_xyz,1)-1,1:3)...
                                =acid_xyz(:,:);
                            counter = counter+size(acid_xyz,1);
                        end
                        if size(acid_fine_xyz,1)>=1
                            resultObject{runCnt}.acid_fine_xyz_accul(counter_f:counter_f + size(acid_fine_xyz,1)-1,1:3)...
                                =acid_fine_xyz(:,:);
                            counter_f = counter_f+size(acid_fine_xyz,1);
                        end
                        if counter_f~=counter
                            fprintf('Length mismatch between pixelated acid positions and the continuous one')
                        end      
                    end
                end
            end
            toc
        elseif exist(fileName_Sep,'file')
            %% Energy separated scan archives
            load(fileName_Sep,'energyScanArchive');            
            if ~exist('energyScanArchive','var')
                disp('Scan Archive exists but not accessible');
                continue
            end
            %%% Find the right incident energy
            tic
            for jj = 1:min(size(energyScanArchive,2),nTrials{runCnt})
                acid_xyz        =   energyScanArchive{jj}.acid_xyz;
                acid_fine_xyz   =   energyScanArchive{jj}.activation_xyz;
                if size(acid_xyz,1)>=1
                    resultObject{runCnt}.acid_xyz_accul(counter:counter + size(acid_xyz,1)-1,1:3)...
                        =acid_xyz(:,:);
                    counter = counter+size(acid_xyz,1);
                end
                if size(acid_fine_xyz,1)>=1
                    resultObject{runCnt}.acid_fine_xyz_accul(counter_f:counter_f + size(acid_fine_xyz,1)-1,1:3)...
                        =acid_fine_xyz(:,:);
                    counter_f = counter_f+size(acid_fine_xyz,1);
                end
                if counter_f~=counter
                    fprintf('Length mismatch between pixelated acid positions and the continuous one')
                end      
            end
            toc           
        else
            for i = 1:nTrials{runCnt}
                %%% Legacy per trial data logging
                try
                    load(strcat('..\..\..\..\',...
                        'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\',...
                        runsToCompare{runCnt},...
                        '\Ein=',incEngy{runCnt},...
                        '_Dose=',doseStr{runCnt},...
                        'epnm2_Ef=15.5_pag-Emin=5_rcnrad=3.00_PAG=0.4_T',num2str(i),'.mat'))
                    if size(acid_xyz,1)>=1
                        resultObject{runCnt}.acid_xyz_accul(counter:counter + size(acid_xyz,1)-1,1:3)...
                            =acid_xyz(:,:);
                        counter = counter+size(acid_xyz,1);
                    end
                    if size(acid_fine_xyz,1)>=1
                        resultObject{runCnt}.acid_fine_xyz_accul(counter_f:counter_f + size(acid_fine_xyz,1)-1,1:3)...
                            =acid_fine_xyz(:,:);
                        counter_f = counter_f+size(acid_fine_xyz,1);
                    end
                    if counter_f~=counter
                        fprintf('Length mismatch between pixelated acid positions and the continuous one')
                    end        
                catch exception
                    disp(exception);
                    fprintf('Problem reading the %d-th file\n', i);
                    break;
                end
            end
        end

        %%% The means
        resultObject{runCnt}.mu_acid =...
            mean(resultObject{runCnt}.acid_xyz_accul);
        resultObject{runCnt}.mu_acid_fine =...
            mean(resultObject{runCnt}.acid_fine_xyz_accul);
        fprintf('The mean of position of acids\n');
        disp(resultObject{runCnt}.mu_acid);

        %%% The spreads
        resultObject{runCnt}.sig_acid =...
            std(resultObject{runCnt}.acid_xyz_accul);
        resultObject{runCnt}.sig_acid_fine =...
            std(resultObject{runCnt}.acid_fine_xyz_accul);
        fprintf('The standard deviation in position of acids\n');
        disp(resultObject{runCnt}.sig_acid);
        fprintf('The standard deviation in position of activation events\n');
        disp(resultObject{runCnt}.sig_acid_fine);    

        %%% Corrected positions
        resultObject{runCnt}.r3=...
            resultObject{runCnt}.acid_xyz_accul...
            -resultObject{runCnt}.mu_acid;
        resultObject{runCnt}.r3_fine=...
            resultObject{runCnt}.acid_fine_xyz_accul...
            -resultObject{runCnt}.mu_acid_fine;

        %%% The radii
        resultObject{runCnt}.rabs =...
            sqrt(sum(resultObject{runCnt}.r3.^2,2));
        resultObject{runCnt}.rabs_fine =...
            sqrt(sum(resultObject{runCnt}.r3_fine.^2,2));

        %%% The event-activation distances
        resultObject{runCnt}.acidEventdist = ...
            sqrt(sum((resultObject{runCnt}.r3_fine-resultObject{runCnt}.r3).^2,2));
    end
end