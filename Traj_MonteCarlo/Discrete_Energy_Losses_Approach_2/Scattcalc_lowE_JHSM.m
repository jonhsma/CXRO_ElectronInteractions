%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What it does: read in a "event trace" which could be incidece (single
%%% event) or a trajectory of an electron. In the latter it determines
%%% which event gives an additional electron and start a new trajectory for
%%% the new electron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eventdata,pagdata,polymdata]=Scattcalc_lowE_JHSM(event,scattdata,eventdata,xyzglobal,pagdata,polymdata,logfile_fid)

global illustration secSpawningTheta;

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
            fprintf(logfile_fid,'WARNING: scatt_Elim is less than optdata_Emin\n');
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
            fprintf(logfile_fid,'.........Scattcalc_lowE: Ese = %.4f eV\n',event{i}.Ese);
            pagimg_pre=pagdata.pagimg;
            %%%% Isotropic Scattering model  
            %[ev2,pagdata,polymdata]=trajcalc3_JHSM_toy(event{i},scattdata,scatt_Elim,xyzglobal,pagdata,polymdata,logfile_fid);
            %%%% Scattering model with proper coordinate transformation
            [ev2,pagdata,polymdata]=trajcalc3_JHSM_angleAdd(event{i},scattdata,scatt_Elim,xyzglobal,pagdata,polymdata,logfile_fid);
  %         fprintf('.........# of PAGs activated = %d\n',sum(pagimg_pre(:))-sum(pagdata.pagimg(:)));
            dbg=1;
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
                    
                    %%%% the emitted secondaries have random angles:
                    % Solid angle compensated distribution
                    cosTheta    =   rand(1)*2-1;
                    ev2{j}.theta_in = acos(cosTheta);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    secSpawningTheta = [secSpawningTheta ev2{j}.theta_in];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ev2{j}.phi_in=2*pi*rand;
                    ev2{j}.scatt_Elim=scatt_Elim;
                    ev2{j}.lowEthr=lowEthr;
                    ev2{j}.lowEimfp=lowEimfp;
                    eventdata{length(eventdata)+1}=ev2{j};
                end
                [eventdata,pagdata,polymdata]=Scattcalc_lowE_JHSM(ev2,scattdata,eventdata,xyzglobal,pagdata,polymdata,logfile_fid);
            end
        end
    end
end