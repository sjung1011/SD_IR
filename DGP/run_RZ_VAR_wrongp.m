%% STATE-DEPENDENT SIMULATION STUDY: Misspecification running file
% SEOJIN JUNG
% This version: June.2025
% Misspecification: 
    % DGPs with a wrong lags p (1 misspec)
    % wrong DGPs (1 misspec)
    % wrong DGPs + wrong lags p (2 misspec)  
% baseline: parameters based on RZ data, short sample, 1e3 MCs 
%% HOUSEKEEPING 

clc
clear all
close all

% Make sure Working directory is C:\Users\seoju\OneDrive - University at Albany - SUNY\SD_25\SD_MAVG_MainResearch
% Add all paths below this folder 
% addpath(genpath(fullfile('..', 'Auxiliary_Functions')));addpath(genpath(fullfile('..','Estimation_Routines')));addpath(genpath('Subroutines'))

rng(1,"twister");
tic;


%% Parallel computing object
num_workers = 4; % number of workers in a parallel pool
poolobj = parpool('local', num_workers);
clear num_workers;

%% SET EXPERIMENT
% check the 'n_MC' in shared.m
exper = 'RZ'; % 'RZ','GHKP'
dgp_type = 'TVAR'; % 'TVAR';'MSVAR'; or 'STVAR' line 173 when misspec is one and two
wrong_dgp = 'TVAR';  % 'TVAR';'MSVAR'; or 'STVAR' % line173 when misspec is two
lag_true = 4; % for dgp
lag_type = 1; % for estimation, 1 and 6; settings.est.nlags
mode_type = 1; % robustness check mode:
               % 1 (nonlinear), 2 (linear), 3 (cumulative IRF),              
sample_size = 'long'; % check line 56, 'long' for 1889:2015 or 'short' for 1947:2015 

%% SETTINGS
%delta_vals = [1,5,10,-1,-5,-10]; % structural shock sizes.

run(fullfile("DGP/Settings/shared.m"));
run(fullfile("DGP/Settings/check_mode.m"));

% Storage folder for results

save_pre = 'wrongS_wrongp_Sta_4mc'; % destination to store the results (test:4, running:120kk)

save_folder = fullfile(save_pre, save_mode_dir,['lag' num2str(lag_type) '_' sample_size]);

%% Call true VAR PARAMETERS
%----------------------------------------------------------------
% load true parameters with lag 4
%----------------------------------------------------------------
% preparation 
T_L = settings.simul.T_L; 
H = settings.est.IRF_hor; % 21
M = settings.simul.n_MC; % share.m 


% extract and store estimated VAR parameters
switch dgp_type
    case 'TVAR'
        file = ['C:\Users\seoju\OneDrive - University at Albany - ' ...
        'SUNY\SD_25\SD_MAVG_MainResearch\Results_Sta_1000mc\nonlinear\lag4_long\RZ_TVAR.mat'];
        load(file, 'DGP_estimate'); 
    case 'MSVAR'
        file = ['C:\Users\seoju\OneDrive - University at Albany - ' ...
        'SUNY\SD_25\SD_MAVG_MainResearch\Results_Sta_1000mc\nonlinear\lag4_long\RZ_MSVAR.mat'];
        load(file, 'DGP_estimate'); 
    case 'STVAR'
        file = ['C:\Users\seoju\OneDrive - University at Albany - ' ...
        'SUNY\SD_25\SD_MAVG_MainResearch\Results_Sta_1000mc\nonlinear\lag4_long\RZ_STVAR.mat'];
        load(file, 'DGP_estimate'); 
end 

%% MONTE CARLO ANALYSIS 
%----------------------------------------------------------------
% Create Placeholders for Results
%----------------------------------------------------------------
% number of estimation methods
settings.est.n_methods = length(settings.est.methods_name); % {'lp','lp_bc','blp','lp_pen','svar','svar_bc','bvar'};
% 
results_irf_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_irf_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_irf_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_irf_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_gic_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_gic_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_gic_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_gic_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_ssr_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_ssr_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_ssr_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_ssr_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_se_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_se_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_se_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_se_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_tss_lps_state0_g = NaN(settings.est.IRF_hor, settings.simul.n_MC); %  H x M, for all lps
results_tss_lps_state1_g = NaN(settings.est.IRF_hor, settings.simul.n_MC); 
results_tss_lps_state0_y = NaN(settings.est.IRF_hor, settings.simul.n_MC); 
results_tss_lps_state1_y = NaN(settings.est.IRF_hor, settings.simul.n_MC); 
results_tss_vars_state0_g = NaN(settings.est.IRF_hor, settings.simul.n_MC); %  H x M, for all vars
results_tss_vars_state1_g = NaN(settings.est.IRF_hor, settings.simul.n_MC); 
results_tss_vars_state0_y = NaN(settings.est.IRF_hor, settings.simul.n_MC); 
results_tss_vars_state1_y = NaN(settings.est.IRF_hor, settings.simul.n_MC);
% CVA
CVA_w_lps_state0_g = NaN(settings.est.IRF_hor,settings.est.n_weights_lp,settings.simul.n_MC); % H x 4 x M
CVA_w_lps_state1_g = NaN(settings.est.IRF_hor,settings.est.n_weights_lp,settings.simul.n_MC);
CVA_w_lps_state0_y = NaN(settings.est.IRF_hor,settings.est.n_weights_lp,settings.simul.n_MC); 
CVA_w_lps_state1_y = NaN(settings.est.IRF_hor,settings.est.n_weights_lp,settings.simul.n_MC);
CVA_w_vars_state0_g = NaN(settings.est.IRF_hor,settings.est.n_weights_var,settings.simul.n_MC); % H x 3 x M
CVA_w_vars_state1_g = NaN(settings.est.IRF_hor,settings.est.n_weights_var,settings.simul.n_MC);
CVA_w_vars_state0_y = NaN(settings.est.IRF_hor,settings.est.n_weights_var,settings.simul.n_MC); 
CVA_w_vars_state1_y = NaN(settings.est.IRF_hor,settings.est.n_weights_var,settings.simul.n_MC);
%
results_hist_state0 = NaN(1,settings.simul.n_MC);
results_hist_state1 = NaN(1,settings.simul.n_MC);

%% --------------------------------------------------------------
% Estimate IRFs, MC loop 
%----------------------------------------------------------------
% preparation 
T_S = settings.simul.T_S; % 250
T_burn =settings.simul.T_burn; %750

%parfor m=1:M % Monte Carlo simulations (parfor location after testing)
parfor m=1:M
    rng(settings.simul.seed(m), 'twister');
    % --------------------------------------------------------------
    % Temporary Storage Folders (due to parfor, inside parfor)
    %---------------------------------------------------------------- 
    temp_irf_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor); % n_method x H 
    temp_irf_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_irf_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_irf_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_gic_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor); 
    temp_gic_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_gic_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_gic_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_ssr_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor); 
    temp_ssr_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_ssr_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_ssr_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_se_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor); 
    temp_se_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_se_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_se_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_tss_lps_state0_g = NaN(settings.est.IRF_hor,1); %  H x 1
    temp_tss_lps_state1_g = NaN(settings.est.IRF_hor,1); 
    temp_tss_lps_state0_y = NaN(settings.est.IRF_hor,1); 
    temp_tss_lps_state1_y = NaN(settings.est.IRF_hor,1); 
    temp_tss_vars_state0_g = NaN(settings.est.IRF_hor,1);
    temp_tss_vars_state1_g = NaN(settings.est.IRF_hor,1); 
    temp_tss_vars_state0_y = NaN(settings.est.IRF_hor,1); 
    temp_tss_vars_state1_y = NaN(settings.est.IRF_hor,1);
    % CVA
    temp_CVA_w_lps_state0_g = NaN(settings.est.IRF_hor,settings.est.n_weights_lp); % H x 4
    temp_CVA_w_lps_state1_g = NaN(settings.est.IRF_hor,settings.est.n_weights_lp);
    temp_CVA_w_lps_state0_y = NaN(settings.est.IRF_hor,settings.est.n_weights_lp); 
    temp_CVA_w_lps_state1_y = NaN(settings.est.IRF_hor,settings.est.n_weights_lp);
    temp_CVA_w_vars_state0_g = NaN(settings.est.IRF_hor,settings.est.n_weights_var); % H x 3
    temp_CVA_w_vars_state1_g = NaN(settings.est.IRF_hor,settings.est.n_weights_var);
    temp_CVA_w_vars_state0_y = NaN(settings.est.IRF_hor,settings.est.n_weights_var); 
    temp_CVA_w_vars_state1_y = NaN(settings.est.IRF_hor,settings.est.n_weights_var);
    %
    temp_hist_state0 = NaN;
    temp_hist_state1 = NaN;

    % PRE-INITIALISE the structural-shock matrices  
    %     (prevents “uninitialised temporary” warnings)
    ei_lp    = [];   ei_bclp  = [];   ei_blp  = [];   ei_pen  = [];
    ei_var   = [];   ei_bcvar = [];   ei_bvar = [];

    % Generate data
    goodDraw = false;
		while ~goodDraw
			try
				[Y_raw, S_raw] = generateData_RZ(T_S,DGP_estimate,lag_true,T_burn,true,wrong_dgp); % Y_raw is w =[x g y]
				goodDraw = true;            % sample is balanced
			catch ME
				if strcmp(ME.identifier,'generateData_RZ:TooFewObs')
					continue                % redraw with same RNG stream
				else
					rethrow(ME)             % any other bug: stop execution
				end
			end
		end
   %
    temp_hist_state0 = sum(S_raw(:) == 0);
    temp_hist_state1 = sum(S_raw(:) == 1);

    % Loop over regimes
    for i_method = 1:settings.est.n_methods
        % estimate IRFs: {'lp','lp_bc','blp','lp_pen','svar','svar_bc','bvar'} % 7
           switch settings.est.methods_name{i_method}
               % lps
               case 'lp' 
                 [IRF, GIC, ei_lp, SSR, TSS, SE]... 
                 = LP_est_stateDep(Y_raw, S_raw, settings, 0); 
                 temp_irf_state0_g(i_method,:,:)=IRF.state0_g;
                 temp_irf_state1_g(i_method,:,:)=IRF.state1_g;
                 temp_irf_state0_y(i_method,:,:)=IRF.state0_y;
                 temp_irf_state1_y(i_method,:,:)=IRF.state1_y;
                 temp_gic_state0_g(i_method,:,:)=GIC.state0_g;
                 temp_gic_state1_g(i_method,:,:)=GIC.state1_g;
                 temp_gic_state0_y(i_method,:,:)=GIC.state0_y;
                 temp_gic_state1_y(i_method,:,:)=GIC.state1_y;
                 temp_ssr_state0_g(i_method,:,:)=SSR.state0_g;
                 temp_ssr_state1_g(i_method,:,:)=SSR.state1_g;
                 temp_ssr_state0_y(i_method,:,:)=SSR.state0_y;
                 temp_ssr_state1_y(i_method,:,:)=SSR.state1_y;
                 temp_se_state0_g(i_method,:,:)=SE.state0_g;
                 temp_se_state1_g(i_method,:,:)=SE.state1_g;
                 temp_se_state0_y(i_method,:,:)=SE.state0_y;
                 temp_se_state1_y(i_method,:,:)=SE.state1_y;
                 temp_tss_lps_state0_g =TSS.state0_g; % only here for lps 
                 temp_tss_lps_state1_g =TSS.state1_g;
                 temp_tss_lps_state0_y =TSS.state0_y;
                 temp_tss_lps_state1_y =TSS.state1_y;

               case 'lp_bc'
                 [IRF, GIC, ei_bclp, SSR, TSS, SE]... 
                 = LP_est_stateDep(Y_raw, S_raw, settings, 1); 
                 temp_irf_state0_g(i_method,:,:)=IRF.state0_g;
                 temp_irf_state1_g(i_method,:,:)=IRF.state1_g;
                 temp_irf_state0_y(i_method,:,:)=IRF.state0_y;
                 temp_irf_state1_y(i_method,:,:)=IRF.state1_y;
                 temp_gic_state0_g(i_method,:,:)=GIC.state0_g;
                 temp_gic_state1_g(i_method,:,:)=GIC.state1_g;
                 temp_gic_state0_y(i_method,:,:)=GIC.state0_y;
                 temp_gic_state1_y(i_method,:,:)=GIC.state1_y;
                 temp_ssr_state0_g(i_method,:,:)=SSR.state0_g;
                 temp_ssr_state1_g(i_method,:,:)=SSR.state1_g;
                 temp_ssr_state0_y(i_method,:,:)=SSR.state0_y;
                 temp_ssr_state1_y(i_method,:,:)=SSR.state1_y;
                 temp_se_state0_g(i_method,:,:)=SE.state0_g;
                 temp_se_state1_g(i_method,:,:)=SE.state1_g;
                 temp_se_state0_y(i_method,:,:)=SE.state0_y;
                 temp_se_state1_y(i_method,:,:)=SE.state1_y;
                 
               case 'blp'
                   [IRF, GIC, ei_blp, SSR, SE]...
                       = BLP_est_stateDep(Y_raw, S_raw, settings);
                 temp_irf_state0_g(i_method,:,:)=IRF.state0_g;
                 temp_irf_state1_g(i_method,:,:)=IRF.state1_g;
                 temp_irf_state0_y(i_method,:,:)=IRF.state0_y;
                 temp_irf_state1_y(i_method,:,:)=IRF.state1_y;
                 temp_gic_state0_g(i_method,:,:)=GIC.state0_g;
                 temp_gic_state1_g(i_method,:,:)=GIC.state1_g;
                 temp_gic_state0_y(i_method,:,:)=GIC.state0_y;
                 temp_gic_state1_y(i_method,:,:)=GIC.state1_y;
                 temp_ssr_state0_g(i_method,:,:)=SSR.state0_g;
                 temp_ssr_state1_g(i_method,:,:)=SSR.state1_g;
                 temp_ssr_state0_y(i_method,:,:)=SSR.state0_y;
                 temp_ssr_state1_y(i_method,:,:)=SSR.state1_y;
                 temp_se_state0_g(i_method,:,:)=SE.state0_g;
                 temp_se_state1_g(i_method,:,:)=SE.state1_g;
                 temp_se_state0_y(i_method,:,:)=SE.state0_y;
                 temp_se_state1_y(i_method,:,:)=SE.state1_y;

               case 'lp_pen'
                   [IRF, GIC, ei_pen, SSR, SE] ...
                    = LP_shrink_est_stateDep(Y_raw, S_raw, settings);
                   temp_irf_state0_g(i_method,:,:)=IRF.state0_g;
                 temp_irf_state1_g(i_method,:,:)=IRF.state1_g;
                 temp_irf_state0_y(i_method,:,:)=IRF.state0_y;
                 temp_irf_state1_y(i_method,:,:)=IRF.state1_y;
                 temp_gic_state0_g(i_method,:,:)=GIC.state0_g;
                 temp_gic_state1_g(i_method,:,:)=GIC.state1_g;
                 temp_gic_state0_y(i_method,:,:)=GIC.state0_y;
                 temp_gic_state1_y(i_method,:,:)=GIC.state1_y;
                 temp_ssr_state0_g(i_method,:,:)=SSR.state0_g;
                 temp_ssr_state1_g(i_method,:,:)=SSR.state1_g;
                 temp_ssr_state0_y(i_method,:,:)=SSR.state0_y;
                 temp_ssr_state1_y(i_method,:,:)=SSR.state1_y;
                 temp_se_state0_g(i_method,:,:)=SE.state0_g;
                 temp_se_state1_g(i_method,:,:)=SE.state1_g;
                 temp_se_state0_y(i_method,:,:)=SE.state0_y;
                 temp_se_state1_y(i_method,:,:)=SE.state1_y;
               % vars
               case 'svar'
                   [IRF, GIC, ei_var, SSR, TSS, SE] = SVAR_est_stateDep(Y_raw, S_raw, settings, 0);
                   temp_irf_state0_g(i_method,:,:)=IRF.state0_g;
                   temp_irf_state1_g(i_method,:,:)=IRF.state1_g;
                   temp_irf_state0_y(i_method,:,:)=IRF.state0_y;
                   temp_irf_state1_y(i_method,:,:)=IRF.state1_y;
                   temp_gic_state0_g(i_method,:,:)=GIC.state0_g;
                   temp_gic_state1_g(i_method,:,:)=GIC.state1_g;
                   temp_gic_state0_y(i_method,:,:)=GIC.state0_y;
                   temp_gic_state1_y(i_method,:,:)=GIC.state1_y;
                   temp_ssr_state0_g(i_method,:,:)=SSR.state0_g;
                   temp_ssr_state1_g(i_method,:,:)=SSR.state1_g;
                   temp_ssr_state0_y(i_method,:,:)=SSR.state0_y;
                   temp_ssr_state1_y(i_method,:,:)=SSR.state1_y;
                   temp_se_state0_g(i_method,:,:)=SE.state0_g;
                   temp_se_state1_g(i_method,:,:)=SE.state1_g;
                   temp_se_state0_y(i_method,:,:)=SE.state0_y;
                   temp_se_state1_y(i_method,:,:)=SE.state1_y;
                   temp_tss_vars_state0_g =TSS.state0_g; % only here for vars 
                   temp_tss_vars_state1_g =TSS.state1_g;
                   temp_tss_vars_state0_y =TSS.state0_y;
                   temp_tss_vars_state1_y =TSS.state1_y;

               case 'svar_bc'
                   [IRF, GIC, ei_bcvar, SSR, TSS, SE] = SVAR_est_stateDep(Y_raw, S_raw, settings, 1);
                   temp_irf_state0_g(i_method,:,:)=IRF.state0_g;
                   temp_irf_state1_g(i_method,:,:)=IRF.state1_g;
                   temp_irf_state0_y(i_method,:,:)=IRF.state0_y;
                   temp_irf_state1_y(i_method,:,:)=IRF.state1_y;
                   temp_gic_state0_g(i_method,:,:)=GIC.state0_g;
                   temp_gic_state1_g(i_method,:,:)=GIC.state1_g;
                   temp_gic_state0_y(i_method,:,:)=GIC.state0_y;
                   temp_gic_state1_y(i_method,:,:)=GIC.state1_y;
                   temp_ssr_state0_g(i_method,:,:)=SSR.state0_g;
                   temp_ssr_state1_g(i_method,:,:)=SSR.state1_g;
                   temp_ssr_state0_y(i_method,:,:)=SSR.state0_y;
                   temp_ssr_state1_y(i_method,:,:)=SSR.state1_y;
                   temp_se_state0_g(i_method,:,:)=SE.state0_g;
                   temp_se_state1_g(i_method,:,:)=SE.state1_g;
                   temp_se_state0_y(i_method,:,:)=SE.state0_y;
                   temp_se_state1_y(i_method,:,:)=SE.state1_y;

               case 'bvar'
                     [IRF, GIC, ei_bvar, SSR, TSS, SE] ...
                    = BVAR_est_stateDep(Y_raw, S_raw, settings);
                   temp_irf_state0_g(i_method,:,:)=IRF.state0_g;
                   temp_irf_state1_g(i_method,:,:)=IRF.state1_g;
                   temp_irf_state0_y(i_method,:,:)=IRF.state0_y;
                   temp_irf_state1_y(i_method,:,:)=IRF.state1_y;
                   temp_gic_state0_g(i_method,:,:)=GIC.state0_g;
                   temp_gic_state1_g(i_method,:,:)=GIC.state1_g;
                   temp_gic_state0_y(i_method,:,:)=GIC.state0_y;
                   temp_gic_state1_y(i_method,:,:)=GIC.state1_y;
                   temp_ssr_state0_g(i_method,:,:)=SSR.state0_g;
                   temp_ssr_state1_g(i_method,:,:)=SSR.state1_g;
                   temp_ssr_state0_y(i_method,:,:)=SSR.state0_y;
                   temp_ssr_state1_y(i_method,:,:)=SSR.state1_y;
                   temp_se_state0_g(i_method,:,:)=SE.state0_g;
                   temp_se_state1_g(i_method,:,:)=SE.state1_g;
                   temp_se_state0_y(i_method,:,:)=SE.state0_y;
                   temp_se_state1_y(i_method,:,:)=SE.state1_y;
           end 
    end % i_method end
    % every ei_* exist compute CVA once
    [CVAw_lp] = CVA_est_lp_stateDep(Y_raw, S_raw, settings,ei_lp,ei_bclp,ei_blp,ei_pen);
    temp_CVA_w_lps_state0_g(:,:) = CVAw_lp.state0_g;
    temp_CVA_w_lps_state1_g(:,:) = CVAw_lp.state1_g;
    temp_CVA_w_lps_state0_y(:,:) = CVAw_lp.state0_y;
    temp_CVA_w_lps_state1_y(:,:) = CVAw_lp.state1_y;

    [CVAw_var] = CVA_est_var_stateDep(Y_raw, S_raw, settings,ei_var,ei_bcvar,ei_bvar);
    temp_CVA_w_vars_state0_g(:,:) = CVAw_var.state0_g;
    temp_CVA_w_vars_state1_g(:,:) = CVAw_var.state1_g;
    temp_CVA_w_vars_state0_y(:,:) = CVAw_var.state0_y;
    temp_CVA_w_vars_state1_y(:,:) = CVAw_var.state1_y;

   % -------------------------------------------
   % Move results to permanat storage in parfor
   % -------------------------------------------
    results_irf_state0_g(:,:,m) = temp_irf_state0_g;
    results_irf_state1_g(:,:,m) = temp_irf_state1_g;
    results_irf_state0_y(:,:,m) = temp_irf_state0_y;
    results_irf_state1_y(:,:,m) = temp_irf_state1_y;
    results_gic_state0_g(:,:,m) = temp_gic_state0_g;
    results_gic_state1_g(:,:,m) = temp_gic_state1_g; 
    results_gic_state0_y(:,:,m) = temp_gic_state0_y; 
    results_gic_state1_y(:,:,m) = temp_gic_state1_y; 
    results_ssr_state0_g(:,:,m) = temp_ssr_state0_g; 
    results_ssr_state1_g(:,:,m) = temp_ssr_state1_g; 
    results_ssr_state0_y(:,:,m) = temp_ssr_state0_y; 
    results_ssr_state1_y(:,:,m) = temp_ssr_state1_y; 
    results_se_state0_g(:,:,m) = temp_se_state0_g; 
    results_se_state1_g(:,:,m) = temp_se_state1_g; 
    results_se_state0_y(:,:,m) = temp_se_state0_y; 
    results_se_state1_y(:,:,m) = temp_se_state1_y; 

    results_tss_lps_state0_g(:,m) = temp_tss_lps_state0_g;
    results_tss_lps_state1_g(:,m) = temp_tss_lps_state1_g; 
    results_tss_lps_state0_y(:,m) = temp_tss_lps_state0_y; 
    results_tss_lps_state1_y(:,m) = temp_tss_lps_state1_y; 
    results_tss_vars_state0_g(:,m) = temp_tss_vars_state0_g;
    results_tss_vars_state1_g(:,m) = temp_tss_vars_state1_g; 
    results_tss_vars_state0_y(:,m) = temp_tss_vars_state0_y; 
    results_tss_vars_state1_y(:,m) = temp_tss_vars_state1_y;
    %
    CVA_w_lps_state0_g(:,:,m) = temp_CVA_w_lps_state0_g;
    CVA_w_lps_state1_g(:,:,m) = temp_CVA_w_lps_state1_g;
    CVA_w_lps_state0_y(:,:,m) = temp_CVA_w_lps_state0_y;
    CVA_w_lps_state1_y(:,:,m) = temp_CVA_w_lps_state1_y;
    CVA_w_vars_state0_g(:,:,m) = temp_CVA_w_vars_state0_g;
    CVA_w_vars_state1_g(:,:,m) = temp_CVA_w_vars_state1_g;
    CVA_w_vars_state0_y(:,:,m) = temp_CVA_w_vars_state0_y;
    CVA_w_vars_state1_y(:,:,m) = temp_CVA_w_vars_state1_y;
    %
    results_hist_state0(m) = temp_hist_state0;
    results_hist_state1(m) = temp_hist_state1;
end 

%% SUMMARIZE RESULTS
% settings.est.IRF_select = 1:21 
 for i_method = 1:settings.est.n_methods
     thisMethod =settings.est.methods_name{i_method};
     % results.irf.(thisMethod).state0_g = squeeze(results_irf_state0_g(i_method,settings.est.IRF_select,:)); % 
     % settings.est.IRF_select =1:21
     results.irf.(thisMethod).state0_g = permute(results_irf_state0_g(i_method,settings.est.IRF_select,:),[2 3 1]); % H x M
     results.irf.(thisMethod).state1_g = permute(results_irf_state1_g(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.irf.(thisMethod).state0_y = permute(results_irf_state0_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.irf.(thisMethod).state1_y = permute(results_irf_state1_y(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.gic.(thisMethod).state0_g = permute(results_gic_state0_g(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.gic.(thisMethod).state1_g = permute(results_gic_state1_g(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.gic.(thisMethod).state0_y = permute(results_gic_state0_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.gic.(thisMethod).state1_y = permute(results_gic_state1_y(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.ssr.(thisMethod).state0_g = permute(results_ssr_state0_g(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.ssr.(thisMethod).state1_g = permute(results_ssr_state1_g(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.ssr.(thisMethod).state0_y = permute(results_ssr_state0_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.ssr.(thisMethod).state1_y = permute(results_ssr_state1_y(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.se.(thisMethod).state0_g = permute(results_se_state0_g(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.se.(thisMethod).state1_g = permute(results_se_state1_g(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.se.(thisMethod).state0_y = permute(results_se_state0_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.se.(thisMethod).state1_y = permute(results_se_state1_y(i_method,settings.est.IRF_select,:),[2 3 1]);

     results.tss.lps.state0_g = squeeze(results_tss_lps_state0_g(settings.est.IRF_select,:));% H x M
     results.tss.lps.state1_g = squeeze(results_tss_lps_state1_g(settings.est.IRF_select,:));
     results.tss.lps.state0_y = squeeze(results_tss_lps_state0_y(settings.est.IRF_select,:));
     results.tss.lps.state1_y = squeeze(results_tss_lps_state1_y(settings.est.IRF_select,:));
     results.tss.vars.state0_g = squeeze(results_tss_vars_state0_g(settings.est.IRF_select,:));     
     results.tss.vars.state1_g = squeeze(results_tss_vars_state1_g(settings.est.IRF_select,:));
     results.tss.vars.state0_y = squeeze(results_tss_vars_state0_y(settings.est.IRF_select,:));
     results.tss.vars.state1_y = squeeze(results_tss_vars_state1_y(settings.est.IRF_select,:));
 end
CVA_w_lps.state0_g = CVA_w_lps_state0_g;
CVA_w_lps.state1_g = CVA_w_lps_state1_g;
CVA_w_lps.state0_y = CVA_w_lps_state0_y;
CVA_w_lps.state1_y = CVA_w_lps_state1_y;
CVA_w_vars.state0_g = CVA_w_vars_state0_g;
CVA_w_vars.state1_g = CVA_w_vars_state1_g;
CVA_w_vars.state0_y = CVA_w_vars_state0_y;
CVA_w_vars.state1_y = CVA_w_vars_state1_y;
%
results.hist.state0 = results_hist_state0;
results.hist.state1 = results_hist_state1;

clear results_* CVA_w_lps_* CVA_w_vars_*
clear i_method m T_burn T_L T_S thisMethod 
%% 
%----------------------------------------------------------------
% Export Results
%----------------------------------------------------------------

mkdir(save_folder);
% e.g., folders: Results_Sta_4mc>nonlinear>lag4 as 'RZ_TVAR_p1.mat' 
save(fullfile(save_folder, strcat(exper, '_', dgp_type,'_', wrong_dgp)), ...
    'DGP_estimate',...
    'settings', ...
    'results', 'CVA_w_lps','CVA_w_vars',...
    'exper','dgp_type','lag_type', 'mode_type','sample_size', ...
    'wrong_dgp','lag_true',...
    '-v7.3'); 

delete(poolobj);
clear save_folder save_pre save_suff save_mode_dir mode_list poolobj

toc;