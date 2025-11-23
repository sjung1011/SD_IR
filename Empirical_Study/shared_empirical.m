%% Empirical settings
% uncertainty shock

%% PREPARATIONS FOR STRUCTURAL ESTIMANDS
% preparation
settings.est.IRF_hor = 37; % (include contemporary)
settings.est.IRF_select = 1:settings.est.IRF_hor;
settings.est.n_state = 2; % expansion, slack

% variable position
var_list   = {'VOLAT', 'STOCK', 'WAGE', 'CPI', 'HOURS', 'EMPM', 'IP','FF'}; % in w
settings.est.resp_vars =[6,7]; % EMPM and IP
settings.est.recurShock_vol = 1; % uncertainty shock
settings.est.normalizeV_vol = 1; % make response to 1 sd volatility shock 
                                 % A shock that moves the S&P index by one unit 
                                 % if it is itslef, a shock that moves itself by one unit
settings.est.recurShock_ff = 8; % monetary shock
settings.est.normalizeV_ff = 8; % make response to 1 pp monetary policy shock

%% ESTIMATION SETTINGS

% set of estimation methods
settings.est.methods_name    = {'svar','bvar','lp','lp_contmp','girf'}; % orthogonal methods
settings.est.nmethods_name = length(settings.est.methods_name);
settings.est.methods_name_linear = {'svar','bvar','lp','lp_contmp'};
settings.est.nmethods_name_linear = length(settings.est.methods_name_linear);

% lag specification
settings.est.n_lag      = 12; % do we use the estimated lag order, or just set it?
settings.est.n_lag_short = 2; % for LPs
% BVAR prior
settings.est.prior.towards_random_walk = 1; % prior shrinking towards random walk? otherwise towards zero
settings.est.bvar_glp                  = 1; % use Giannone, Lenza & Primiceri (2015) BVAR procedure?
                                            % otherwise use basic BVAR with default MN prior (see settings below)
settings.est.prior.tight_overall       = 0.04;
settings.est.prior.tight_nonown_lag    = 0.25;
settings.est.prior.decay_power         = 2;
settings.est.prior.tight_exogenous     = 1e5;

% BVAR posterior draws: raise in empirical

settings.est.posterior_ndraw = 5e3; % number of posterior draws (if set to 0, only use posterior mean of VAR coefficient to compute posterior mean of IRF)

% LP smoothing

settings.est.lambdaRange   = [0.001:0.005:0.021, 0.05:0.1:1.05, 2:1:19, 20:20:100, 200:200:2000]; % cross validation grid, scaled up by T
settings.est.irfLimitOrder = 2; % shrink towards polynomial of that order
settings.est.CV_folds      = 5; % Number of folds used for cross validation

% Optimal weights
%settings.est.wmethods_name = {'EQ', 'BMA', 'JMA'}; % weigths methods  
%settings.est.lpestm_name = {'lp'};
settings.est.n_weights_lp = 2; % for LP approaches
settings.est.n_weights_var = 2; % for VAR approaches

% bootstrap, other draws
settings.est.spec_btstrp = 'wild'; % options: 'standard', 'wild'
settings.est.n_btstrp = 5e3; 

% girf
settings.est.n_girf_paths =5e3;
settings.est.delta_mode = 'onesd'; % 'unit','onesd','twosd'