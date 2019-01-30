%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function isolate all photoemission events and store it in an array
% to save space. Photoemission simulations are gigantic and the resultant
% energyScanArchive is usually huge. To avoid saving data in matlab 7.3
% files, which are usually  huge, this function is written
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [escapeEvents,incidentEvents,preceedingEvents] = extractEscapeEvents(energyScanArchive)
    % find out the ceiling for memory usage
    escapeEvents(1) = energyScanArchive{1}.incidences{1}{1};
    incidentEvents(1) = escapeEvents(1);
    preceedingEvents(1) = escapeEvents(1);
    counter = 1;
    loopCounter = 0;
    terminationCounter = 0;
    for ii = 1:size(energyScanArchive,2)
        for jj =  1:size(energyScanArchive{ii}.incidences,2)
            for kk =  1:size(energyScanArchive{ii}.incidences{jj},2)
                loopCounter = loopCounter+1;
                if strcmp(energyScanArchive{ii}.incidences{jj}{kk}.act,...
                        'escape')
                    escapeEvents(counter) = energyScanArchive{ii}.incidences{jj}{kk};
                    incidentEvents(counter) = energyScanArchive{ii}.incidences{jj}{1};
                    if kk > 1 &&...
                            ~strcmp(energyScanArchive{ii}.incidences{jj}{kk-1}.act,...
                            'escape')&&...
                            ~strcmp(energyScanArchive{ii}.incidences{jj}{kk-1}.act,...
                            'StoneWall')                        
                        preceedingEvents(counter) = energyScanArchive{ii}.incidences{jj}{kk-1};
                    else
                        preceedingEvents(counter) = energyScanArchive{ii}.incidences{jj}{kk};
                        preceedingEvents(counter).scattType = 'NewTraj';
                    end
                    
                    if kk<length(energyScanArchive{ii}.incidences{jj})&&...
                            sum(energyScanArchive{ii}.incidences{jj}{kk}.xyz ==...
                           energyScanArchive{ii}.incidences{jj}{kk+1}.xyz_init)>0
                       disp('Warning! trajectory continues after escape')
                    end
                    counter = counter+1;
                    terminationCounter = terminationCounter + 1;
                elseif strcmp(energyScanArchive{ii}.incidences{jj}{kk}.act,...
                        'StoneWall')
                    terminationCounter = terminationCounter + 1;
                end
            end
        end
    end   
    disp(loopCounter)
    disp(terminationCounter)
    sprintf('Total number of electrons simulated : %d\n',terminationCounter);
end