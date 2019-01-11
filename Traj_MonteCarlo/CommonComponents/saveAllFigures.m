%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function saves all existing graphs as .pngs to a given folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = saveAllFigures(outputFolder)
    y = 0;
    if ~exist(strcat(outputFolder, '\graphics'),'file')
           mkdir(strcat(outputFolder, '\graphics'));
    end
    graphicFolderName = strcat(outputFolder, '\graphics');   
    figList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(figList)
      figHandle = figList(iFig);
      figName   = get(figHandle, 'Name');
      %savefig(figHandle, fullfile(graphicFolderName, figName, '.fig'));
      if ~isempty(num2str(figHandle.Number))>0
        saveas(figHandle,fullfile(graphicFolderName,...
            strcat('figure',num2str(figHandle.Number),...
            figName,'.png')));
      end
    end
    % if the function returns one it means success
    y = 1;
end