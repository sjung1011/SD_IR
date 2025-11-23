%% DGP settings

%% PREPARATIONS FOR STRUCTURAL ESTIMANDS

% preparation
settings.est.IRF_hor = 21; % (include contemporary)
settings.est.IRF_select = 1:settings.est.IRF_hor;
settings.est.n_state = 2; % expansion, slack
% settings.est.useLowerTriang = false; % false: C is not lower triangle

% variable position
settings.est.shock = 1; % index of the shock 
settings.est.resp_vars = [2,3]; % indices of response variables of interest 

% number of Monte Carlo draws
settings.simul.n_MC    = 1000; % number of Monte Carlo reps 1e3 (robustness 5e3)
settings.simul.seed    = (1:settings.simul.n_MC)*10 + randi([0,9],1,settings.simul.n_MC); % random seed for each Monte Carlo
%settings.simul.seed    = 06132025;

% sample settings
settings.simul.T_L     = 2000; % time periods for averaging population IR (potential outcomes)

settings.simul.T_S      = 1500; % time periods for prediction outcomes in each simulation
settings.simul.T_burn = 500; % burn-in

%% ESTIMATION SETTINGS

% set of estimation methods

settings.est.methods_name    = {'lp','lp_bc','blp','lp_pen','svar','svar_bc','bvar'};
%{'lp','lp_bc','blp','lp_pen','svar','svar_bc','bvar'}; % 7

% lag specification

settings.est.n_lag      = lag_type; % do we use the estimated lag order, or just set it?

% BVAR prior

settings.est.prior.towards_random_walk = 1; % prior shrinking towards random walk? otherwise towards zero
settings.est.bvar_glp                  = 1; % use Giannone, Lenza & Primiceri (2015) BVAR procedure?
                                            % otherwise use basic BVAR with default MN prior (see settings below)
settings.est.prior.tight_overall       = 0.04;
settings.est.prior.tight_nonown_lag    = 0.25;
settings.est.prior.decay_power         = 2;
settings.est.prior.tight_exogenous     = 1e5;

% BVAR posterior draws

settings.est.posterior_ndraw = 100; % number of posterior draws (if set to 0, only use posterior mean of VAR coefficient to compute posterior mean of IRF)

% LP smoothing

settings.est.lambdaRange   = [0.001:0.005:0.021, 0.05:0.1:1.05, 2:1:19, 20:20:100, 200:200:2000]; % cross validation grid, scaled up by T
settings.est.irfLimitOrder = 2; % shrink towards polynomial of that order
settings.est.CV_folds      = 5; % Number of folds used for cross validation

% Optimal weights
settings.est.case_name = {'case1', 'case2'}; % case1: 1:4 LALP, case2: four LP estimation methods
settings.est.wmethods_name = {'EQ', 'BMA', 'JMA'}; % weigths methods  
settings.est.lpestm_name = {'lp','blp','lp_corrbias','lp_penalize'}; %case2
settings.est.n_weights_lp = 4; % for LP approaches
settings.est.n_weights_var = 3; % for VAR approaches

% GIRF
settings.est.n_girf_paths = 1000; % over 300
settings.est.delta_mode = 'onesd'; % 'unit', 'onesd' or 'twosd'