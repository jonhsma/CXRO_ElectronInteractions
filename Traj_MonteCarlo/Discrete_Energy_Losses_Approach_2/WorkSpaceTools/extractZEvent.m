%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function isolate all events that crosses the specified z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [passageEvents,incidentEvents] = extractZEvent(energyScanArchive,z,direction)
    passageEvents(1) = energyScanArchive{1}.incidences{1}{1};
    incidentEvents(1) = passageEvents(1);
    counter = 1;
    loopCounter = 0;
    terminationCounter = 0;
    
    % Make sure direction is always +- 1
    if direction~=0
        sgn = direction./abs(direction);
    else
        sgn = 0;
    end
    
    for ii = 1:size(energyScanArchive,2)
        for jj =  1:size(energyScanArchive{ii}.incidences,2)
            for kk =  1:size(energyScanArchive{ii}.incidences{jj},2)
                
                loopCounter = loopCounter+1;
                
                %initial and final z of the event
                z_i = energyScanArchive{ii}.incidences{jj}{kk}.xyz_init(3);
                z_f = energyScanArchive{ii}.incidences{jj}{kk}.xyz(3);
                
                % The condition
                if sgn ==0
                    condition = (z_i < z && z_f > z) ||(-z_i < -z && -z_f > -z);
                else
                    condition = sgn*z_i < sgn*z && sgn*z_f > sgn*z;
                end
                
                if condition
                    passageEvents(counter) = energyScanArchive{ii}.incidences{jj}{kk};
                    incidentEvents(counter) = energyScanArchive{ii}.incidences{jj}{1};
                    counter = counter+1;
                end
            end
        end
    end   
    disp(loopCounter)
    sprintf('Total number of electrons simulated : %d\n',terminationCounter);
    if counter == 1
        % counter is always grater than the number of valid events by 1
        passageEvents = [];
        incidentEvents = [];
    end
end