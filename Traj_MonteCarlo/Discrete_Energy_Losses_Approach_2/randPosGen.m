%%% This function generates the positions of a certain species randomly. It
%%% is mathematically guaranteed that the distribution of number of that
%%% species per unit volume is a Poisson Distribution. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% dimension   (1,:) = [lowerlimit upperlimit] in x direction
%%%             (2,:) = [lowerlimit upperlimit] in y direction (so on)
%%% nExpected   the expectation value of the TOTAL number of the
%%%             species IN THE ENTIRE VOLUME. The actually total number
%%%             will be poissrnd(nExpected);
%%%             If number is negative then the total number will be
%%%             nExpected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function positions = randPosGen(limits,nExpected)
    %% The total number of the species is a Poisson number as well
    if nExpected > 0 
        nTotal = poissrnd(nExpected);
    else
        nTotal = -nExpected;
    end    
    %% Some rearrangements
    gridSize = limits(:,2)-limits(:,1);
    %% Generating the posisions
    positions = rand([3 nTotal]).*gridSize+limits(:,1);    
end