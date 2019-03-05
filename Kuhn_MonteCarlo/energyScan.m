% Thsi script runs the monte carlo for the given energies 
function results = energyScan(energyArray,nSegments,scattdata,varargin)

    % output parameters for length analysis
    Lcoef       =   zeros([2 length(energyArray)]);
    Lcoef_se    =   zeros([2 length(energyArray)]);
    % output parameters for length analysis
    Lcoef_ms       =   zeros([2 length(energyArray)]);
    Lcoef_se_ms    =   zeros([2 length(energyArray)]);
    % The expectation values
    Lmean       =   zeros([nSegments length(energyArray)]);
    Lse         =   Lmean;
    % The expectation values of R^2
    Lms         =   Lmean;
    Lsse        =   Lmean;
    
    if nargin == 4 % A new follower is selected
        follower = varargin{1};
    else
        follower = @TrajectoryFollowerConstEnergy;
    end
    
    parfor ii = 1:length(energyArray)
        % Initial configuration
        initConfig = struct();
        initConfig.energy = -1;
        initConfig.xyz = [0 0 0];
        initConfig.theta_in = 0; % initial condition for theta in z
        initConfig.phi_in = 0; % initial condition for the polar angle
        initConfig.nSteps = nSegments;  
        initConfig.energy = energyArray(ii);
        
        % The trajectory
        currentTraj = follower(initConfig,...
            scattdata,...
            12/120*6.02*1e23,...
            @arrayAllocator);
        
        % Length analysis
        [sum_L,~] = trajLengthAna(currentTraj);        
        range = 30:400;
        % Fit the distance as a function of step
        error = (sum_L.ini.std)./sqrt(sum_L.ini.N+1);
        lm = fitlm(log(range),log(sum_L.ini.mean(range)),...
            'Weights',sum_L.ini.mean(range)./error(range));
        coef = lm.Coefficients.Estimate;
        % Record the coeffcients
        Lcoef(:,ii)     = coef;
        Lcoef_se(:,ii)  = lm.Coefficients.SE; 
        % Fit the distance squared as a function of step
        error = (sum_L.ini.sstd)./sqrt(sum_L.ini.N+1);
        lm = fitlm(range,sum_L.ini.ms(range),...
            'Weights',1./error(range),...
            'Intercept',false);
        coef = lm.Coefficients.Estimate;
        % Record the coeffcients
        Lcoef_ms(:,ii)     = coef;
        Lcoef_se_ms(:,ii)  = lm.Coefficients.SE; 
        

        % Record the expectation values
        Lmean(:,ii)     = sum_L.ini.mean;
        Lse(:,ii)       = sum_L.ini.std;  
        
        Lms(:,ii)       = sum_L.ini.ms;
        Lsse(:,ii)      = sum_L.ini.sstd;
    end
    
    % Save the results
    results.Lcoef       = Lcoef;
    results.Lcoef_se    = Lcoef_se;
    results.Lcoef_ms    = Lcoef_ms;
    results.Lcoef_se_ms = Lcoef_se_ms;
    % The step - end2end distributions
    results.Lmean       = Lmean;
    results.Lsd         = Lse;
    results.Lms         = Lms;
    results.Lsdse       = Lsse;
end