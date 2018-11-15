%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function displays the trajectory of a trial in the scanArchive
%%% database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = displayTrajectories_Trial(database,E_count,trIdx)
    %%  Extract the relavant trial
    trialData       =   database{E_count}{trIdx};    
    dosingSequence  =   trialData.dosingSequence;
    acid_xyz        =   trialData.acid_xyz;
    acid_act_e_xyz  =   trialData.activation_xyz;

    %%  Configure the graphics
    FIGURE_TRAJ_PER_TRIAL   =   7201; % Temporary
    tdHandle = figure(FIGURE_TRAJ_PER_TRIAL);
    tdHandle.Name = 'TrajMonitor';
    grid on
    hold off
    
    database = 1;
    
    %% Creating the graph
    if size(dosingSequence,2) <= 1
        %% Fewer than 20 trajectories
        lengenEntries = cell(size(dosingSequence,2)+2);
        for inIterator = 1:size(dosingSequence,2)
            qX = [];
            qY = [];
            qZ = [];
            qU = [];
            qV = [];
            qW = [];
            nEvents     =   length(...
                trialData.incidences{inIterator});
            for eventIterator = 1:nEvents
                this = trialData.incidences{inIterator}{eventIterator};
                qX = [qX this.xyz_init(1)];
                qY = [qY this.xyz_init(2)];
                qZ = [qZ this.xyz_init(3)];
                qU = [qU this.xyz(1)];
                qV = [qV this.xyz(2)];
                qW = [qW this.xyz(3)];
            end
            quiver3(qX,qY,qZ,qU-qX,qV-qY,qW-qZ,0,'-','LineWidth',2,...
                'Color',hsv2rgb(...
                [inIterator/size(dosingSequence,2),...
                1,1]));
            hold on
            legendEntries{inIterator} = strcat('Electron ',num2str(inIterator));
        end
        if length(acid_xyz) >=1
            quiver3(acid_act_e_xyz(:,1),acid_act_e_xyz(:,2),acid_act_e_xyz(:,3),...
                acid_xyz(:,1) - acid_act_e_xyz(:,1),...
                acid_xyz(:,2) - acid_act_e_xyz(:,2),...
                acid_xyz(:,3) - acid_act_e_xyz(:,3),...
                0,'--');
            scatter3(acid_xyz(:,1),acid_xyz(:,2),acid_xyz(:,3),'o')            
            legendEntries{size(dosingSequence,2)+1} = 'Acid Activation';
            legendEntries{size(dosingSequence,2)+2} = 'Acids';
        end
    else
        %% More than 20 trajectories
        qX = [];
        qY = [];
        qZ = [];
        qU = [];
        qV = [];
        qW = [];
        for inIterator = 1:size(trialData.incidences,2)
            nEvents     =   length(...
                trialData.incidences{inIterator});
            for eventIterator = 1:nEvents
                this = trialData.incidences{inIterator}{eventIterator};
                qX = [qX this.xyz_init(1)];
                qY = [qY this.xyz_init(2)];
                qZ = [qZ this.xyz_init(3)];
                qU = [qU this.xyz(1)];
                qV = [qV this.xyz(2)];
                qW = [qW this.xyz(3)];
            end
        end
        quiver3(qX,qY,qZ,qU-qX,qV-qY,qW-qZ,0,'.-','LineWidth',1);
        legendEntries   =   {'Trajectories'};
        hold on            
        if length(acid_xyz) >=1
            quiver3(acid_act_e_xyz(:,1),acid_act_e_xyz(:,2),acid_act_e_xyz(:,3),...
                acid_xyz(:,1) - acid_act_e_xyz(:,1),...
                acid_xyz(:,2) - acid_act_e_xyz(:,2),...
                acid_xyz(:,3) - acid_act_e_xyz(:,3),...
                0,'--');
            scatter3(acid_xyz(:,1),acid_xyz(:,2),acid_xyz(:,3),'o') 
            legendEntries = [legendEntries {'Activations'},{'Acids Activated'}];
        end
    end
    legend(legendEntries);
    title('Trajectories and activation events in the last iteration')
    daspect([1 1 1])
    drawnow;
    y = tdHandle
end