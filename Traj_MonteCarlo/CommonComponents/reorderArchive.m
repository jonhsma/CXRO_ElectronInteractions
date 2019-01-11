%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reorders an energyScanArchive into a strcuture containing
% no cell array. This is a measure to increase execution/analysis speed
% and save sapce and more importnatly loading time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some notes
%--------------------------------------------------------------------------
% 1 Don't use the sum of nSE as the total number of SE. A secondary only
% materializes (and be initialized in the wrapper) if it's energy is above
% cutoff. So the sum of nSE is always greater than the number of secondary
% electron trajectories
%--------------------------------------------------------------------------
function database = reorderArchive(energyScanArchive)
    
    nTrials = length(energyScanArchive);
    eventCounter = 0;
    incidenceCounter = 0;
    
    % Created for immediate access
    secondaryCounter = 0;
    terminationCounter = 0;
    
    trialsInfoArray      =   [];
    incidenceInfoArray =    [];    
    
    % Extract the length of the events array for initializing the array
    % storing all events
    nTotalEvents = 0;
    for ii = 1:nTrials
        nIncidences = length(energyScanArchive{ii}.incidences);
        for jj = 1: nIncidences
            nEvents = length(energyScanArchive{ii}.incidences{jj});
            nTotalEvents = nTotalEvents + nEvents;
        end
    end
    eventsArray(nEvents)    =   energyScanArchive{1}.incidences{1}{1};
    
    for ii = 1:nTrials
        nIncidences = length(energyScanArchive{ii}.incidences);
        for jj = 1: nIncidences
            incidenceCounter = incidenceCounter + 1;
            nEvents = length(energyScanArchive{ii}.incidences{jj});
            incidenceInfo.secondaryStartingEventIndex  =   [];
            incidenceInfo.terminationEventIndex        =   [];
            for kk = 1:nEvents
                    eventCounter = eventCounter+1;
                    % if a new trial starts
                    if jj == 1 && kk == 1
                        trialInfo.startingEventIndex = eventCounter;
                        trialInfo.startingIncidenceIndex = incidenceCounter;
                    end
                    % if a new incidence starts
                    if kk == 1
                        incidenceInfo.startingEventIndex = eventCounter;
                    end
                    eventsArray(eventCounter) = energyScanArchive{ii}.incidences{jj}{kk};
                    
                    %if kk >1 && sum((energyScanArchive{ii}.incidences{jj}{kk}.xyz_init...
                     %       -energyScanArchive{ii}.incidences{jj}{kk-1}.xyz).^2) > 1e-8
                    if kk >1 && energyScanArchive{ii}.incidences{jj}{kk}.Ein ~=...
                            energyScanArchive{ii}.incidences{jj}{kk-1}.Eout
                    
                        % Register the indices of these events
                        incidenceInfo.terminationEventIndex =...
                            [incidenceInfo.terminationEventIndex eventCounter-1];
                        incidenceInfo.secondaryStartingEventIndex  =...
                            [incidenceInfo.secondaryStartingEventIndex eventCounter];
                        % The through trial counter for easy asccess
                        secondaryCounter    = secondaryCounter + 1;
                        terminationCounter  = terminationCounter + 1;                        
                    end
            end
            %include the last termination event
            incidenceInfo.terminationEventIndex =...
                    [incidenceInfo.terminationEventIndex eventCounter];
            terminationCounter  = terminationCounter + 1;
            
            incidenceInfo.endingEventIndex = eventCounter;
            incidenceInfo.trialNumber = ii;
            incidenceInfoArray = [incidenceInfoArray incidenceInfo];
        end
        trialInfo.endingEventIndex  =   eventCounter;
        trialInfo.endingIncidenceIndex  =   incidenceCounter;
        trialInfo.dosingSequence   =   energyScanArchive{ii}.dosingSequence;
        trialInfo.acid_xyz         =   energyScanArchive{ii}.acid_xyz;
        trialInfo.activation_xyz   =   energyScanArchive{ii}.activation_xyz;
        trialsInfoArray = [trialsInfoArray trialInfo];
        %disp('Trial completed');
    end
    
    database.eventsArray        =   eventsArray;
    database.incidenceInfoArray =   incidenceInfoArray;
    database.trialsInfoArray    =   trialsInfoArray;
    database.nTermination       =   terminationCounter;
    database.nSecondary         =   secondaryCounter;
end