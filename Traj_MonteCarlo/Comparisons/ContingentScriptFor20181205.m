%% Contigency Loading Script For 20181205_Ref
diary on
diary('20181206_for20181205_Ref')
list = {'scanArchive_20',...
    'scanArchive_30',...
    'scanArchive_45',...
    'scanArchive_50',...
    'scanArchive_65',...
    'scanArchive_80'};


resultObject = cell(1,length(list));

for ii = 1:length(list)
    counter = 1;
    counter_f = 1;
    load(strcat('..\..\..\..\',...
    'JonathanCodeIO_CXRO\ElectronInteractions\LEEMRes\',...
    '20181205_Ref_SemiColFixed','\',list{ii},'.mat'));
    eval(strcat('scanArchive =',list{ii},';')) 
    for jj = 1:min(size(scanArchive{1},2),1100)
        acid_xyz        =   scanArchive{1}{jj}.acid_xyz;
        acid_fine_xyz   =   scanArchive{1}{jj}.activation_xyz;
        if size(acid_xyz,1)>=1
            resultObject{ii}.acid_xyz_accul(counter:counter + size(acid_xyz,1)-1,1:3)...
                =acid_xyz(:,:);
            counter = counter+size(acid_xyz,1);
        end
        if size(acid_fine_xyz,1)>=1
            resultObject{ii}.acid_fine_xyz_accul(counter_f:counter_f + size(acid_fine_xyz,1)-1,1:3)...
                =acid_fine_xyz(:,:);
            counter_f = counter_f+size(acid_fine_xyz,1);
        end
        if counter_f~=counter
            fprintf('Length mismatch between pixelated acid positions and the continuous one')
        end      
    end
    
    %%% The means
    resultObject{ii}.mu_acid =...
        mean(resultObject{ii}.acid_xyz_accul);
    resultObject{ii}.mu_acid_fine =...
        mean(resultObject{ii}.acid_fine_xyz_accul);
    fprintf(strcat('The mean of position of acids for ',list{ii},'\n'));
    disp(resultObject{ii}.mu_acid);

    %%% The spreads
    resultObject{ii}.sig_acid =...
        std(resultObject{ii}.acid_xyz_accul);
    resultObject{ii}.sig_acid_fine =...
        std(resultObject{ii}.acid_fine_xyz_accul);
    fprintf('The standard deviation in position of acids\n');
    disp(resultObject{ii}.sig_acid);
    fprintf('The standard deviation in position of activation events\n');
    disp(resultObject{ii}.sig_acid_fine); 

end
