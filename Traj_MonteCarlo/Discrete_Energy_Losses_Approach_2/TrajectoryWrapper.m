%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What it does: read in a "event trace" which could be incidece (single
%%% event) or a trajectory of an electron. In the latter it determines
%%% which event gives an additional electron and start a new trajectory for
%%% the new electron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eventdata,pagdata,polymdata]=TrajectoryWrapper(event,scattdata,eventdata,xyzglobal,pagdata,polymdata,logfile_fid)

global illustration secSpawningTheta echoConfig demoEcho;

    if isfield(scattdata,'followerHandle') &&...
            isa(scattdata.followerHandle,'function_handle')
        trajectoryFollower  =   scattdata.followerHandle;
    else
        trajectoryFollower  =   @TrajectoryFollower;
    end
    
   
        

    plotscatt=0;
    % if isfield(scattdata.optical,'inel_dcsdata')
    %     scatt_Elim=min(scattdata.optical.E);
    % else
    %     scatt_Elim=3; % how low of energy to trak the electrons at?
    % end
    % scatt_Elim=20; % over ride the above if-else

    scatt_Elim=event{1}.scatt_Elim;
    lowEthr=event{1}.lowEthr;
    lowEimfp=event{1}.lowEimfp;

    if isfield(scattdata.optical,'inel_dcsdata')
        optdata_Emin=min(scattdata.optical.E);
        if scatt_Elim<optdata_Emin
            %fprintf(logfile_fid,'WARNING: scatt_Elim is less than optdata_Emin\n');
        end
    end

    %%%% global co-ordinates [needed for modeling forward scattering]
    % xyzglobal.x=event.x;
    % xyzglobal.y=event.y;
    % xyzglobal.z=event.z;

    for i = 1:length(event)
        if illustration
            fprintf('Iteration %d in ScattCalc_lowE\n',i);
        end
        if event{i}.Ese>scatt_Elim
            %fprintf(logfile_fid,'.........Scattcalc_lowE: Ese = %.4f eV\n',event{i}.Ese);
            posPAG_init     =   pagdata.posPAG;
            %% Trajectory Follower
            [ev2,pagdata,polymdata]=trajectoryFollower(event{i},scattdata,scatt_Elim,xyzglobal,pagdata,polymdata,logfile_fid);
            %% Analyse the trajectory and generate secondary trajectories
            if demoEcho
                fprintf('# of PAGs activated in the last execution of the follower = %d\n',...
                   sum(isnan(pagdata.posPAG(1,:)))-sum(isnan(posPAG_init(1,:))));
            end
            % A dummy event array for secondary generation
            % so that nothing in ev2 is altered. ev2 can thus serve as a
            % accurate record of the trajectory
            evInc = ev2;            
            if ~isempty(ev2)
                for j=1:length(ev2)
                    xyzglobal.x=[xyzglobal.x ev2{j}.xyz(1)];
                    xyzglobal.y=[xyzglobal.y ev2{j}.xyz(2)];
                    xyzglobal.z=[xyzglobal.z ev2{j}.xyz(3)];
                    if plotscatt~=0
                        figure(fig1);
                        hold on;
        %                 plot(ev2{j}.x,ev2{j}.y,'o');drawnow;
                        plot3(ev2{j}.xyz(1),ev2{j}.xyz(2),ev2{j}.xyz(3),'o');drawnow;
                        xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
                    end
                    %%% All events go to eventdata
                    eventdata{length(eventdata)+1}=ev2{j};
                    
                    %% Secondary seeding
                    %%%% the emitted secondaries have random angles:
                    % Solid angle compensated distribution
                    %.theta_in (and .phi_in) means the theta of the photoelectron generated
                    cosTheta    =   rand(1)*2-1;
                    evInc{j}.theta_in = acos(cosTheta); 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    secSpawningTheta = [secSpawningTheta evInc{j}.theta_in];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    evInc{j}.phi_in=2*pi*rand;
                    evInc{j}.scatt_Elim=scatt_Elim;
                    evInc{j}.lowEthr=lowEthr;
                    evInc{j}.lowEimfp=lowEimfp;                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% eventdata, pagdata and polymdata are "threading
                %%% variables" They are updated everytime the follower is
                %%% run. The input of one iteration (of i) is the outcome
                %%% of the previous one. If no scattering happens these
                %%% variables threads through the functions without any
                %%% change.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [eventdata,pagdata,polymdata]=TrajectoryWrapper(evInc,scattdata,eventdata,xyzglobal,pagdata,polymdata,logfile_fid);
            end
        end
    end
end