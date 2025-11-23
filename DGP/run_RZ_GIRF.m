%% STATE-DEPENDENT SIMULATION STUDY: GIRF running FILE
% SEOJIN JUNG
% This version: July.2025
% baseline: parameters based on RZ data, short sample, 1e3 MCs 
% DGPs based on a transition process  
% GIRF_VAR(p)s use VAR(p)'s GIC  
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
dgp_type = {'TVAR','MSVAR','STVAR'}; 
lag_type = 4; % No. of lags to impose in estimation
mode_type = 1; % robustness check mode:
               % 1 (nonlinear), 2 (linear), 3 (cumulative IRF),              
sample_size = 'long'; % check line 56, 'long' for 1889:2015 or 'short' for 1947:2015 
data_trans = 'stationary'; % 'level', 'stationary'
%% SETTINGS
% delta_vals = 1; % [1,5,10,-1,-5,-10]; % structural shock sizes.

run(fullfile("DGP/Settings/shared.m"));
run(fullfile("DGP/Settings/check_mode.m"));

% Storage folder for results

save_pre = 'GIRF_Sta_1000mc'; % destination to store the results (test:4, running:120kk)
save_folder = fullfile(save_pre, save_mode_dir,['lag' num2str(lag_type) '_' sample_size]);

%% DGP: VAR PARAMETERS should be inside parfor
%----------------------------------------------------------------
% load data and estimate parameters
%----------------------------------------------------------------
% preparation 
T_L = settings.simul.T_L; 
H = settings.est.IRF_hor; % 21
M = settings.simul.n_MC; % share.m 

% extract and store estimated VAR parameters
% TVAR
file = ['C:\Users\seoju\OneDrive - University at Albany - ' ...
        'SUNY\SD_25\SD_MAVG_MainResearch\Results_Sta_1000mc\nonlinear\lag4_long\RZ_TVAR.mat'];
load(file, 'DGP_estimate');
DGP_estimate_TVAR = DGP_estimate;
% MSVAR
file = ['C:\Users\seoju\OneDrive - University at Albany - ' ...
        'SUNY\SD_25\SD_MAVG_MainResearch\Results_Sta_1000mc\nonlinear\lag4_long\RZ_MSVAR.mat'];
load(file, 'DGP_estimate');
DGP_estimate_MSVAR = DGP_estimate;
% STVAR
file = ['C:\Users\seoju\OneDrive - University at Albany - ' ...
        'SUNY\SD_25\SD_MAVG_MainResearch\Results_Sta_1000mc\nonlinear\lag4_long\RZ_STVAR.mat'];
load(file, 'DGP_estimate');
DGP_estimate_STVAR = DGP_estimate;

%% Parallel computing object: 
num_workers = 4; % number of workers in a parallel pool
poolobj = parpool('local', num_workers);
clear num_workers;
%% MONTE CARLO ANALYSIS 
%----------------------------------------------------------------
% Create Placeholders for Results
%----------------------------------------------------------------
% number of estimation methods
settings.est.n_methods = length(dgp_type); 
results_irf_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_irf_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_irf_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_irf_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_bandlo_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_bandlo_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_bandlo_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_bandlo_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_bandhi_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); % n_method x H x M
results_bandhi_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_bandhi_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 
results_bandhi_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor, settings.simul.n_MC); 

%{
CVA_w_vars_state0_g = NaN(settings.est.IRF_hor,settings.est.n_methods,settings.simul.n_MC); % H x 3 x M
CVA_w_vars_state1_g = NaN(settings.est.IRF_hor,settings.est.n_methods,settings.simul.n_MC);
CVA_w_vars_state0_y = NaN(settings.est.IRF_hor,settings.est.n_methods,settings.simul.n_MC); 
CVA_w_vars_state1_y = NaN(settings.est.IRF_hor,settings.est.n_methods,settings.simul.n_MC);
%}
results_hist_state0 = NaN(settings.est.n_methods,settings.simul.n_MC);
results_hist_state1 = NaN(settings.est.n_methods,settings.simul.n_MC);
%% --------------------------------------------------------------
% Estimate IRFs, MC loop 
%----------------------------------------------------------------
% preparation 
T_S = settings.simul.T_S; % 400
T_burn =settings.simul.T_burn; %600

parfor m=1:M % Monte Carlo simulations (parfor location after testing)
    rng(settings.simul.seed(m), 'twister');
    % --------------------------------------------------------------
    % Temporary Storage Folders (due to parfor, inside parfor)
    %---------------------------------------------------------------- 
    temp_irf_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor); % n_method x H 
    temp_irf_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_irf_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_irf_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_bandlo_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor); 
    temp_bandlo_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_bandlo_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_bandlo_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_bandhi_state0_g = NaN(settings.est.n_methods, settings.est.IRF_hor); 
    temp_bandhi_state1_g = NaN(settings.est.n_methods, settings.est.IRF_hor);
    temp_bandhi_state0_y = NaN(settings.est.n_methods, settings.est.IRF_hor);  
    temp_bandhi_state1_y = NaN(settings.est.n_methods, settings.est.IRF_hor);
   
    %{ 
    temp_CVA_w_vars_state0_g = NaN(settings.est.IRF_hor,settings.est.n_methods); % H x 3
    temp_CVA_w_vars_state1_g = NaN(settings.est.IRF_hor,settings.est.n_methods);
    temp_CVA_w_vars_state0_y = NaN(settings.est.IRF_hor,settings.est.n_methods); 
    temp_CVA_w_vars_state1_y = NaN(settings.est.IRF_hor,settings.est.n_methods);
    %}
    temp_hist_state0 = NaN(settings.est.n_methods,1); % 3 x 1
    temp_hist_state1 = NaN(settings.est.n_methods,1);
    % PRE-INITIALISE the structural-shock matrices  
    %{     
    ei_tvar   = [];   ei_msvar = [];   ei_stvar = [];
    %}
    for i_method = 1:length(dgp_type)
        switch dgp_type{i_method}
            case 'TVAR'
                % Generate data
	            goodDraw = false;
		        while ~goodDraw
			            try
				[Y_raw, S_raw] = generateData_RZ(T_S,DGP_estimate_TVAR,lag_type,T_burn,true,'TVAR'); % Y_raw is w =[x g y]
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
                temp_hist_state0(i_method) = sum(S_raw(:) == 0);
                temp_hist_state1(i_method) = sum(S_raw(:) == 1);

                % fcn 
                [IRF, band_lo, band_hi, ~] = GIRF_VAR_est_stateDep(Y_raw, S_raw, settings);
                temp_irf_state0_g(i_method,:,:)=IRF.state0_g;
                temp_irf_state1_g(i_method,:,:)=IRF.state1_g;
                temp_irf_state0_y(i_method,:,:)=IRF.state0_y;
                temp_irf_state1_y(i_method,:,:)=IRF.state1_y;
                temp_bandlo_state0_g(i_method,:,:)=band_lo.state0_g;
                temp_bandlo_state1_g(i_method,:,:)=band_lo.state1_g;
                temp_bandlo_state0_y(i_method,:,:)=band_lo.state0_y;
                temp_bandlo_state1_y(i_method,:,:)=band_lo.state1_y;
                temp_bandhi_state0_g(i_method,:,:)=band_hi.state0_g;
                temp_bandhi_state1_g(i_method,:,:)=band_hi.state1_g;
                temp_bandhi_state0_y(i_method,:,:)=band_hi.state0_y;
                temp_bandhi_state1_y(i_method,:,:)=band_hi.state1_y;
             
            case 'MSVAR'
                % Generate data
	            goodDraw = false;
		        while ~goodDraw
			            try
				[Y_raw, S_raw] = generateData_RZ(T_S,DGP_estimate_MSVAR,lag_type,T_burn,true,'MSVAR'); % Y_raw is w =[x g y]
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
                temp_hist_state0(i_method) = sum(S_raw(:) == 0);
                temp_hist_state1(i_method) = sum(S_raw(:) == 1);

                % fcn 
                [IRF, band_lo, band_hi, ~] = GIRF_VAR_est_stateDep(Y_raw, S_raw, settings);
                temp_irf_state0_g(i_method,:,:)=IRF.state0_g;
                temp_irf_state1_g(i_method,:,:)=IRF.state1_g;
                temp_irf_state0_y(i_method,:,:)=IRF.state0_y;
                temp_irf_state1_y(i_method,:,:)=IRF.state1_y;
                temp_bandlo_state0_g(i_method,:,:)=band_lo.state0_g;
                temp_bandlo_state1_g(i_method,:,:)=band_lo.state1_g;
                temp_bandlo_state0_y(i_method,:,:)=band_lo.state0_y;
                temp_bandlo_state1_y(i_method,:,:)=band_lo.state1_y;
                temp_bandhi_state0_g(i_method,:,:)=band_hi.state0_g;
                temp_bandhi_state1_g(i_method,:,:)=band_hi.state1_g;
                temp_bandhi_state0_y(i_method,:,:)=band_hi.state0_y;
                temp_bandhi_state1_y(i_method,:,:)=band_hi.state1_y;

            case 'STVAR'
                % Generate data
	            goodDraw = false;
		        while ~goodDraw
			            try
				[Y_raw, S_raw] = generateData_RZ(T_S,DGP_estimate_STVAR,lag_type,T_burn,true,'STVAR'); % Y_raw is w =[x g y]
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
                temp_hist_state0(i_method) = sum(S_raw(:) == 0);
                temp_hist_state1(i_method) = sum(S_raw(:) == 1);

                % fcn 
                [IRF, band_lo, band_hi, ~] = GIRF_VAR_est_stateDep(Y_raw, S_raw, settings);
                temp_irf_state0_g(i_method,:,:)=IRF.state0_g;
                temp_irf_state1_g(i_method,:,:)=IRF.state1_g;
                temp_irf_state0_y(i_method,:,:)=IRF.state0_y;
                temp_irf_state1_y(i_method,:,:)=IRF.state1_y;
                temp_bandlo_state0_g(i_method,:,:)=band_lo.state0_g;
                temp_bandlo_state1_g(i_method,:,:)=band_lo.state1_g;
                temp_bandlo_state0_y(i_method,:,:)=band_lo.state0_y;
                temp_bandlo_state1_y(i_method,:,:)=band_lo.state1_y;
                temp_bandhi_state0_g(i_method,:,:)=band_hi.state0_g;
                temp_bandhi_state1_g(i_method,:,:)=band_hi.state1_g;
                temp_bandhi_state0_y(i_method,:,:)=band_hi.state0_y;
                temp_bandhi_state1_y(i_method,:,:)=band_hi.state1_y;
        end
    end
 
    %{ 
    % technically impossible to use this function as Y, S will differ
    depends on transition process. 
    [CVAw_var] = CVA_est_var_stateDep(Y_raw, S_raw, settings,ei_tvar,ei_msvar,ei_stvar);
    temp_CVA_w_vars_state0_g(:,:) = CVAw_var.state0_g;
    temp_CVA_w_vars_state1_g(:,:) = CVAw_var.state1_g;
    temp_CVA_w_vars_state0_y(:,:) = CVAw_var.state0_y;
    temp_CVA_w_vars_state1_y(:,:) = CVAw_var.state1_y;
    %}

   % -------------------------------------------
   % Move results to permanat storage in parfor
   % -------------------------------------------
    results_irf_state0_g(:,:,m) = temp_irf_state0_g;
    results_irf_state1_g(:,:,m) = temp_irf_state1_g;
    results_irf_state0_y(:,:,m) = temp_irf_state0_y;
    results_irf_state1_y(:,:,m) = temp_irf_state1_y;
    results_bandlo_state0_g(:,:,m) = temp_bandlo_state0_g;
    results_bandlo_state1_g(:,:,m) = temp_bandlo_state1_g; 
    results_bandlo_state0_y(:,:,m) = temp_bandlo_state0_y; 
    results_bandlo_state1_y(:,:,m) = temp_bandlo_state1_y; 
    results_bandhi_state0_g(:,:,m) = temp_bandhi_state0_g; 
    results_bandhi_state1_g(:,:,m) = temp_bandhi_state1_g; 
    results_bandhi_state0_y(:,:,m) = temp_bandhi_state0_y; 
    results_bandhi_state1_y(:,:,m) = temp_bandhi_state1_y; 
    %{
    CVA_w_vars_state0_g(:,:,m) = temp_CVA_w_vars_state0_g;
    CVA_w_vars_state1_g(:,:,m) = temp_CVA_w_vars_state1_g;
    CVA_w_vars_state0_y(:,:,m) = temp_CVA_w_vars_state0_y;
    CVA_w_vars_state1_y(:,:,m) = temp_CVA_w_vars_state1_y;
    %}
    results_hist_state0(:,m) = temp_hist_state0;
    results_hist_state1(:,m) = temp_hist_state1;
end 


%% SUMMARIZE RESULTS
% settings.est.IRF_select = 1:21 
 for i_method = 1:settings.est.n_methods
     thisMethod = dgp_type{i_method};
     % results.irf.(thisMethod).state0_g = squeeze(results_irf_state0_g(i_method,settings.est.IRF_select,:)); % 
     % settings.est.IRF_select =1:21
     results.irf.(thisMethod).state0_g = permute(results_irf_state0_g(i_method,settings.est.IRF_select,:),[2 3 1]); % H x M
     results.irf.(thisMethod).state1_g = permute(results_irf_state1_g(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.irf.(thisMethod).state0_y = permute(results_irf_state0_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.irf.(thisMethod).state1_y = permute(results_irf_state1_y(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.bandlo.(thisMethod).state0_g = permute(results_bandlo_state0_g(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.bandlo.(thisMethod).state1_g = permute(results_bandlo_state1_g(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.bandlo.(thisMethod).state0_y = permute(results_bandlo_state0_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.bandlo.(thisMethod).state1_y = permute(results_bandlo_state1_y(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.bandhi.(thisMethod).state0_g = permute(results_bandhi_state0_g(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.bandhi.(thisMethod).state1_g = permute(results_bandhi_state1_g(i_method,settings.est.IRF_select,:),[2 3 1]);
     results.bandhi.(thisMethod).state0_y = permute(results_bandhi_state0_y(i_method,settings.est.IRF_select,:),[2 3 1]); 
     results.bandhi.(thisMethod).state1_y = permute(results_bandhi_state1_y(i_method,settings.est.IRF_select,:),[2 3 1]);
     % 
     results.hist.(thisMethod).state0 = results_hist_state0(i_method,:);
     results.hist.(thisMethod).state1 = results_hist_state1(i_method,:);
 end
%{
CVA_w_vars.state0_g = CVA_w_vars_state0_g;
CVA_w_vars.state1_g = CVA_w_vars_state1_g;
CVA_w_vars.state0_y = CVA_w_vars_state0_y;
CVA_w_vars.state1_y = CVA_w_vars_state1_y;
%}

clear results_* CVA_w_vars_*
clear H i_method m M T_burn T_L T_S thisMethod 
%% 
%----------------------------------------------------------------
% Export Results
%----------------------------------------------------------------

mkdir(save_folder);
% e.g., folders: Results_Sta_4mc>nonlinear>lag4 as 'RZ_TVAR_p1.mat' 
save(fullfile(save_folder, strcat(exper, '_allDGPs')), ...
    'settings', ...
    'results',...
    'exper','dgp_type','lag_type', 'mode_type','sample_size', 'data_trans',...
    '-v7.3'); %'CVA_w_vars',

delete(poolobj);
clear save_folder save_pre save_suff save_mode_dir mode_list poolobj

toc;
