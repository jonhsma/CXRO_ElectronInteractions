function outCat = weightedCategoricalRandGen(cats,weights)
    if length(cats) ~= length(weights)
        outCat  = 'Warning';
    end
    cumWeight   =   cumsum(weights);
    cumWeight   =   cumWeight/cumWeight(end);
    
    dice        =   rand(1);
    
    outCat      =   cats{find(dice<cumWeight,1)};
end