%% STATE-DEPENDENT SIMULATION STUDY: linear process for comparison
% SEOJIN JUNG
% This version: June.2025
% Parameters based on RZ data, short sample, 1e3 MCs 
% DGPs based on a transition process  
%% HOUSEKEEPING 

clc
clear all
close all

% Make sure Working directory is C:\Users\seoju\OneDrive - University at Albany - SUNY\SD_25\SD_MAVG_MainResearch
% Add all paths below this folder 
% addpath(genpath(fullfile('..', 'Auxiliary_Functions')));addpath(genpath(fullfile('..','Estimation_Routines')));addpath(genpath('Subroutines'))
rng(1, 'twister');
tic;

%% SET EXPERIMENT
% check the 'n_MC' in shared.m
exper = 'RZ'; % 'RZ','GHKP'
dgp_type = 'LSVAR'; % No state
lag_type = 4; % No. of lags to impose in estimation
mode_type = 2; % robustness check mode:
               % 1 (nonlinear), 2 (linear), 3 (cumulative IRF),              
sample_size = 'long'; % check line 56, 'long' for 1889:2015 or 'short' for 1947:2015 

%% SETTINGS
delta_vals = [1,5,10,-1,-5,-10]; % structural shock sizes.

run(fullfile("DGP/Settings/shared.m"));
run(fullfile("DGP/Settings/check_mode.m"));

% Storage folder for results

save_pre = 'Results_Sta_1000mc'; % destination to store the results (test:4, running:120kk)

%save_suff = num2str(lag_type);
%save_folder = fullfile(save_pre, save_mode_dir, strcat('lag', save_suff));
save_folder = fullfile(save_pre, save_mode_dir,['lag' num2str(lag_type) '_' sample_size]);
%% DGP: VAR PARAMETERS
%----------------------------------------------------------------
% load data and estimate parameters of linear process
%----------------------------------------------------------------
% preparation 
T_L = settings.simul.T_L; 

% extract and store estimated VAR parameters
switch exper   
    case 'RZ' % RZ (2018) for fiscal narrative instrument shock experiment
        DGP_estimate = estimateParams_RZ_linear('RZDAT.xlsx','short','LSVAR');
end

%% MONTE CARLO SIMULATIONS for CARs
%----------------------------------------------------------------
% Generate Data based on DGP_estimate and compute CARs 
% the true model is LS VAR (4)
%----------------------------------------------------------------
CAR_g = zeros(settings.est.IRF_hor, length(delta_vals)); % 21x6
CAR_y = zeros(settings.est.IRF_hor, length(delta_vals));
H = settings.est.IRF_hor; % 21
M = settings.simul.n_MC; % share.m 

%rng(settings.simul.seed, 'twister');
for m = 1: M
    rng(settings.simul.seed(m), 'twister');
    % Generate new dataset for each Monte Carlo iteration
    [synthetic_w] = generateData_RZ_linear(T_L,DGP_estimate,lag_type,0,false);

    % Temporary storage for impulse responses for this dataset for each regime i
    CAR_g_i = zeros(H, length(delta_vals));
    CAR_y_i = zeros(H, length(delta_vals));
  
    R_i = size(synthetic_w,1); % Number of obs
        % Loop over shock sizes
        for d = 1:length(delta_vals)
            delta = delta_vals(d);
            num_skipped = 0; % Reset skipped count for each delta
            
            % Temporary storage for responses across histories for the current delta
            CAR_g_temp = zeros(H, R_i);
            CAR_y_temp = zeros(H, R_i);
            valid_index = 0; % Track valid cases separately

            for t = 1:R_i % Iterate over time
                
                if t > 3
                   % Use the previous four lag information for histories to generate potential output, 4x3
                   Y_hist4 = synthetic_w(t-3:t , :);
                else
                    % If t ≤ 3, throw away this data to compute CARs
                    num_skipped = num_skipped + 1;
                    continue
                end
                
                % Increment valid cases counter
                valid_index = valid_index + 1;

                % Path 1: With structural shock
                Y_with_shock = simulatePath_linear(Y_hist4, DGP_estimate, delta);

                % Path 2: Without structural shock
                Y_without_shock = simulatePath_linear(Y_hist4, DGP_estimate, 0);

                % Compute the difference
                CAR_g_temp(:, valid_index) = Y_with_shock(:, 2) - Y_without_shock(:, 2); % Government spending response
                CAR_y_temp(:, valid_index) = Y_with_shock(:, 3) - Y_without_shock(:, 3); % GDP response
            end
            
            % Average across valid histories for regime i
            valid_cases = R_i - num_skipped;
            if valid_cases > 0
                CAR_g_i(:, d) = sum(CAR_g_temp(:, 1:valid_cases), 2) / valid_cases;
                CAR_y_i(:, d) = sum(CAR_y_temp(:, 1:valid_cases), 2) / valid_cases;
            end
        end 
end
    % Accumulate impulse responses across datasets (for each MC)
    CAR_g = CAR_g + CAR_g_i;
    CAR_y = CAR_y + CAR_y_i;


% Average across all datasets
CAR_g = CAR_g / M;
CAR_y = CAR_y / M;

clear CAR_g_i CAR_g_temp CAR_y_i CAR_y_temp d delta num_skipped synthetic_S ... 
synthetic_w t valid_cases valid_index Y_hist4 Y_without_shock Y_with_shock

%% Parallel computing object

num_workers = 4; % number of workers in a parallel pool
poolobj = parpool('local', num_workers);
clear num_workers;

%% MONTE CARLO ANALYSIS 
% Misspecification: estimates with linear process 
% Compare it with CARs of non-linear process, not the CARs of linear process 
% ----------------------------------------------------------------
% Create Placeholders for Results
%----------------------------------------------------------------
% number of estimation methods
settings.est.n_methods = length(settings.est.methods_name); % {'lp','lp_bc','blp','lp_pen','svar','svar_bc','bvar'};
% 
results_irf_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_irf_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_gic_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_gic_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_ssr_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_ssr_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_se_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_se_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_tss_lps_g = NaN(settings.est.IRF_hor, settings.simul.n_MC); %  H x M, for all lps
results_tss_lps_y = NaN(settings.est.IRF_hor, settings.simul.n_MC); 
results_tss_vars_g = NaN(settings.est.IRF_hor, settings.simul.n_MC); %  H x M, for all vars
results_tss_vars_y = NaN(settings.est.IRF_hor, settings.simul.n_MC); 
% CVA
CVA_w_lps_g = NaN(settings.est.IRF_hor,settings.est.n_weights_lp,settings.simul.n_MC); % H x 4 x M
CVA_w_lps_y = NaN(settings.est.IRF_hor,settings.est.n_weights_lp,settings.simul.n_MC); 
CVA_w_vars_g = NaN(settings.est.IRF_hor,settings.est.n_weights_var,settings.simul.n_MC); % H x 3 x M
CVA_w_vars_y = NaN(settings.est.IRF_hor,settings.est.n_weights_var,settings.simul.n_MC); 

%% --------------------------------------------------------------
% Estimate IRFs, MC loop 
%----------------------------------------------------------------
% preparation 
T_S = settings.simul.T_S; % 250
T_burn =settings.simul.T_burn; %750


parfor m=1:M % Monte Carlo simulations (parfor location after testing)
    rng(settings.simul.seed(m), 'twister');
    % --------------------------------------------------------------
    % Temporary Storage Folders (due to parfor, inside parfor)
    %---------------------------------------------------------------- 
    temp_irf_g = NaN(settings.est.n_methods, settings.est.IRF_hor); % n_method x H 
    temp_irf_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_gic_g = NaN(settings.est.n_methods, settings.est.IRF_hor); 
    temp_gic_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_ssr_g = NaN(settings.est.n_methods, settings.est.IRF_hor); 
    temp_ssr_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_se_g = NaN(settings.est.n_methods, settings.est.IRF_hor); 
    temp_se_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_tss_lps_g = NaN(settings.est.IRF_hor,1); %  H x 1
    temp_tss_lps_y = NaN(settings.est.IRF_hor,1); 
    temp_tss_vars_g = NaN(settings.est.IRF_hor,1);
    temp_tss_vars_y = NaN(settings.est.IRF_hor,1); 
    % CVA
    temp_CVA_w_lps_g = NaN(settings.est.IRF_hor,settings.est.n_weights_lp); % H x 4
    temp_CVA_w_lps_y = NaN(settings.est.IRF_hor,settings.est.n_weights_lp); 
    temp_CVA_w_vars_g = NaN(settings.est.IRF_hor,settings.est.n_weights_var); % H x 3
    temp_CVA_w_vars_y = NaN(settings.est.IRF_hor,settings.est.n_weights_var); 

    % settings for estimate fcns
    respV_g = settings.est.resp_vars(1);
    respV_y = settings.est.resp_vars(2);
    recurShock = settings.est.shock; 
    normalizeV = recurShock;

    
    % PRE-INITIALISE the structural-shock matrices  
    ei_lp_g    = [];   ei_bclp_g  = [];   ei_blp_g  = [];   ei_pen_g  = [];
    ei_var_g   = [];   ei_bcvar_g = [];   ei_bvar_g = [];
    ei_lp_y    = [];   ei_bclp_y  = [];   ei_blp_y  = [];   ei_pen_y  = [];
    ei_var_y   = [];   ei_bcvar_y = [];   ei_bvar_y = [];
    
    % Generate data
    [Y_raw] = generateData_RZ_linear(T_S,DGP_estimate,lag_type,T_burn,true); % Y_raw is w =[x g y]
    
    % Loop over regimes
    for i_method = 1:settings.est.n_methods
        % estimate IRFs: {'lp','lp_bc','blp','lp_pen','svar','svar_bc','bvar'} % 7
           switch settings.est.methods_name{i_method}
               % lps
               case 'lp' 
                 [IRF_g, GIC_g, ei_lp_g, SSR_g, TSS_g, SE_g]...
                    = LP_est(Y_raw, settings, respV_g,recurShock, normalizeV,0);  
                 temp_irf_g(i_method,:,:)=IRF_g;
                 temp_gic_g(i_method,:,:)=GIC_g;
                 temp_ssr_g(i_method,:,:)=SSR_g;
                 temp_se_g(i_method,:,:)=SE_g;
                 temp_tss_lps_g =TSS_g; % only here for lps 

                 [IRF_y, GIC_y, ei_lp_y, SSR_y, TSS_y, SE_y]...
                    = LP_est(Y_raw, settings,respV_y,recurShock, normalizeV,0);
                 temp_irf_y(i_method,:,:)=IRF_y;
                 temp_gic_y(i_method,:,:)=GIC_y;
                 temp_ssr_y(i_method,:,:)=SSR_y;
                 temp_se_y(i_method,:,:)=SE_y;
                 temp_tss_lps_y =TSS_y;
               
               case 'lp_bc'
                   [IRF_g, GIC_g, ei_bclp_g, SSR_g, TSS_g, SE_g]...
                    = LP_est(Y_raw, settings,respV_g,recurShock, normalizeV,1);  
                 temp_irf_g(i_method,:,:)=IRF_g;
                 temp_gic_g(i_method,:,:)=GIC_g;
                 temp_ssr_g(i_method,:,:)=SSR_g;
                 temp_se_g(i_method,:,:)=SE_g;
                
                 [IRF_y, GIC_y, ei_bclp_y, SSR_y, TSS_y, SE_y]...
                    = LP_est(Y_raw, settings,respV_y,recurShock, normalizeV,1);
                 temp_irf_y(i_method,:,:)=IRF_y;
                 temp_gic_y(i_method,:,:)=GIC_y;
                 temp_ssr_y(i_method,:,:)=SSR_y;
                 temp_se_y(i_method,:,:)=SE_y;
                 
               case 'blp'
                   [IRF_g, GIC_g, ei_blp_g, SSR_g, SE_g]... 
                   = BLP_est(Y_raw, settings, respV_g, normalizeV);
                 temp_irf_g(i_method,:,:)=IRF_g;
                 temp_gic_g(i_method,:,:)=GIC_g;
                 temp_ssr_g(i_method,:,:)=SSR_g;
                 temp_se_g(i_method,:,:)=SE_g;

                   [IRF_y, GIC_y, ei_blp_y, SSR_y, SE_y]... 
                   = BLP_est(Y_raw, settings, respV_y, normalizeV);
                 temp_irf_y(i_method,:,:)=IRF_y;
                 temp_gic_y(i_method,:,:)=GIC_y;
                 temp_ssr_y(i_method,:,:)=SSR_y;
                 temp_se_y(i_method,:,:)=SE_y;

               case 'lp_pen'
                   [IRF_g, GIC_g, ei_pen_g, SSR_g, SE_g, ~] ...
                       = LP_shrink_est(Y_raw,settings ,respV_g,recurShock,normalizeV);
                 temp_irf_g(i_method,:,:)=IRF_g;
                 temp_gic_g(i_method,:,:)=GIC_g;
                 temp_ssr_g(i_method,:,:)=SSR_g;
                 temp_se_g(i_method,:,:)=SE_g;

                   [IRF_y, GIC_y, ei_pen_y, SSR_y, SE_y, ~]...
                       = LP_shrink_est(Y_raw,settings ,respV_y,recurShock,normalizeV);
                 temp_irf_y(i_method,:,:)=IRF_y;
                 temp_gic_y(i_method,:,:)=GIC_y;
                 temp_ssr_y(i_method,:,:)=SSR_y;
                 temp_se_y(i_method,:,:)=SE_y;

               % vars
               case 'svar'
                   [IRF_g,GIC_g, ei_var_g, SSR_g, TSS_g, SE_g] ...
                       = SVAR_est(Y_raw,settings,respV_g,recurShock,normalizeV,0);
                 temp_irf_g(i_method,:,:)=IRF_g;
                 temp_gic_g(i_method,:,:)=GIC_g;
                 temp_ssr_g(i_method,:,:)=SSR_g;
                 temp_se_g(i_method,:,:)=SE_g;
                 temp_tss_vars_g =TSS_g; % only here for vars

                 [IRF_y,GIC_y, ei_var_y, SSR_y, TSS_y, SE_y] ...
                       = SVAR_est(Y_raw,settings,respV_y,recurShock,normalizeV,0);
                 temp_irf_y(i_method,:,:)=IRF_y;
                 temp_gic_y(i_method,:,:)=GIC_y;
                 temp_ssr_y(i_method,:,:)=SSR_y;
                 temp_se_y(i_method,:,:)=SE_y;
                 temp_tss_vars_y =TSS_y; 

               case 'svar_bc'
                   [IRF_g,GIC_g, ei_bcvar_g, SSR_g, TSS_g, SE_g] ...
                       = SVAR_est(Y_raw,settings,respV_g,recurShock,normalizeV,1);
                 temp_irf_g(i_method,:,:)=IRF_g;
                 temp_gic_g(i_method,:,:)=GIC_g;
                 temp_ssr_g(i_method,:,:)=SSR_g;
                 temp_se_g(i_method,:,:)=SE_g;

                 [IRF_y,GIC_y, ei_bcvar_y, SSR_y, TSS_y, SE_y] ...
                       = SVAR_est(Y_raw,settings,respV_y,recurShock,normalizeV,1);
                 temp_irf_y(i_method,:,:)=IRF_y;
                 temp_gic_y(i_method,:,:)=GIC_y;
                 temp_ssr_y(i_method,:,:)=SSR_y;
                 temp_se_y(i_method,:,:)=SE_y;

               case 'bvar'
                   [IRF_g,GIC_g,ei_bvar_g,SSR_g,SE_g]...
                       = BVAR_est(Y_raw,settings, respV_g,recurShock,normalizeV);
                 temp_irf_g(i_method,:,:)=IRF_g;
                 temp_gic_g(i_method,:,:)=GIC_g;
                 temp_ssr_g(i_method,:,:)=SSR_g;
                 temp_se_g(i_method,:,:)=SE_g;

                 [IRF_y,GIC_y,ei_bvar_y,SSR_y,SE_y]...
                       = BVAR_est(Y_raw,settings, respV_y,recurShock,normalizeV);
                 temp_irf_y(i_method,:,:)=IRF_y;
                 temp_gic_y(i_method,:,:)=GIC_y;
                 temp_ssr_y(i_method,:,:)=SSR_y;
                 temp_se_y(i_method,:,:)=SE_y;
           end 
    end % i_method end
    % every ei_* exist compute CVA once
    [CVAw_lp_g] ...
        = CVA_est_lp(Y_raw,settings,ei_lp_g,ei_bclp_g,ei_blp_g,ei_pen_g, respV_g,recurShock);
    temp_CVA_w_lps_g(:,:) = CVAw_lp_g';
    [CVAw_lp_y] ...
        = CVA_est_lp(Y_raw,settings,ei_lp_y,ei_bclp_y,ei_blp_y,ei_pen_y, respV_y,recurShock);
    temp_CVA_w_lps_y(:,:) = CVAw_lp_y';
    [CVAw_var_g] ...
        = CVA_est_var(Y_raw, settings,ei_var_g,ei_bcvar_g,ei_bvar_g,respV_g);
    temp_CVA_w_vars_g(:,:) = CVAw_var_g';
    [CVAw_var_y] ...
        = CVA_est_var(Y_raw, settings,ei_var_y,ei_bcvar_y,ei_bvar_y,respV_y);
    temp_CVA_w_vars_y(:,:) = CVAw_var_y';

   % -------------------------------------------
   % Move results to permanat storage in parfor
   % -------------------------------------------
    results_irf_g(:,:,m) = temp_irf_g;
    results_irf_y(:,:,m) = temp_irf_y;
    results_gic_g(:,:,m) = temp_gic_g;
    results_gic_y(:,:,m) = temp_gic_y;  
    results_ssr_g(:,:,m) = temp_ssr_g; 
    results_ssr_y(:,:,m) = temp_ssr_y; 
    results_se_g(:,:,m) = temp_se_g; 
    results_se_y(:,:,m) = temp_se_y; 
    results_tss_lps_g(:,m) = temp_tss_lps_g;
    results_tss_lps_y(:,m) = temp_tss_lps_y; 
    results_tss_vars_g(:,m) = temp_tss_vars_g;
    results_tss_vars_y(:,m) = temp_tss_vars_y; 
    %
    CVA_w_lps_g(:,:,m) = temp_CVA_w_lps_g;
    CVA_w_lps_y(:,:,m) = temp_CVA_w_lps_y;
    CVA_w_vars_g(:,:,m) = temp_CVA_w_vars_g;
    CVA_w_vars_y(:,:,m) = temp_CVA_w_vars_y;
end 

%% SUMMARIZE RESULTS
% settings.est.IRF_select = 1:21 
 for i_method = 1:settings.est.n_methods
     thisMethod =settings.est.methods_name{i_method};
     % results.irf.(thisMethod).state0_g = squeeze(results_irf_state0_g(i_method,settings.est.IRF_select,:)); % 
     % settings.est.IRF_select =1:21
     results.irf.(thisMethod).g = permute(results_irf_g(i_method,settings.est.IRF_select,:),[2 3 1]); % H x M
     results.irf.(thisMethod).y = permute(results_irf_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.gic.(thisMethod).g = permute(results_gic_g(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.gic.(thisMethod).y = permute(results_gic_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.ssr.(thisMethod).g = permute(results_ssr_g(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.ssr.(thisMethod).y = permute(results_ssr_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.se.(thisMethod).g = permute(results_se_g(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.se.(thisMethod).y = permute(results_se_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     
     results.tss.lps.g = squeeze(results_tss_lps_g(settings.est.IRF_select,:));% H x M
     results.tss.lps.y = squeeze(results_tss_lps_y(settings.est.IRF_select,:));
     results.tss.vars.g = squeeze(results_tss_vars_g(settings.est.IRF_select,:));     
     results.tss.vars.y = squeeze(results_tss_vars_y(settings.est.IRF_select,:));
 end
CVA_w_lps.g = CVA_w_lps_g;
CVA_w_lps.y = CVA_w_lps_y;
CVA_w_vars.g = CVA_w_vars_g;
CVA_w_vars.y = CVA_w_vars_y;

clear results_* CVA_w_lps_* CVA_w_vars_*
clear H i i_method m M r R_i T_burn T_L T_S thisMethod 
%% 
%----------------------------------------------------------------
% Export Results
%----------------------------------------------------------------

mkdir(save_folder);
% e.g., folders: Results_Sta_4mc>nonlinear>lag4 as 'RZ_TVAR_p1.mat' 
save(fullfile(save_folder, strcat(exper, '_', dgp_type)), ...
    'CAR_g', 'CAR_y','DGP_estimate','delta_vals',...
    'settings', ...
    'results', 'CVA_w_lps','CVA_w_vars',...
    'exper','dgp_type','lag_type', 'mode_type','sample_size', ...
     ...
    '-v7.3'); 

delete(poolobj);
clear save_folder save_pre save_suff save_mode_dir mode_list poolobj

toc;
