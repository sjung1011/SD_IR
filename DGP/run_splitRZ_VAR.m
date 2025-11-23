%% STATE-DEPENDENT SIMULATION STUDY I: MAIN FILE
% SEOJIN JUNG
% This version: June.2025
% split data first then generate data for each regime. 
    % Not a big different in CARs bw original and this version

%% HOUSEKEEPING 

clc
clear all
close all

% Make sure Working directory is C:\Users\seoju\OneDrive - University at Albany - SUNY\SD_25\SD_MAVG_MainResearch
% Add all paths below this folder 
% addpath(genpath(fullfile('..', 'Auxiliary_Functions')));addpath(genpath(fullfile('..','Estimation_Routines')));addpath(genpath('Subroutines'))

%% SET EXPERIMENT
% check the 'n_MC' in shared.m
exper = 'GHKP'; % 'RZ','GHKP'
dgp_type = 'TVAR'; % 'TVAR';'MAVAR';'MSVAR'; or 'STVAR'
lag_type = 4; % No. of lags to impose in estimation
mode_type = 1; % robustness check mode:
               % 1 (nonlinear), 2 (linear), 3 (cumulative IRF),              
sample_size = 'short'; % check line 56, 'long' for 1889:2015 or 'short' for 1947:2015
%% SETTINGS
delta_vals = [1,5,10,-1,-5,-10]; % structural shock sizes.

run(fullfile("DGP/Settings/shared.m"));
run(fullfile("DGP/Settings/check_mode.m"));

% Storage folder for results

save_pre = 'Split_Sta_4mc'; % destination to store the results (test:4, running:120kk)

save_suff = num2str(lag_type);

save_folder = fullfile(save_pre, save_mode_dir, strcat('lag', save_suff));

%% DGP: VAR PARAMETERS
%----------------------------------------------------------------
% load data and estimate parameters
%----------------------------------------------------------------
% preparation 
T_L = settings.simul.T_L; 

% extract and store estimated VAR parameters
switch exper   
    case 'RZ' % parameterize based on RZ (2018) for fiscal narrative instrument shock experiment
        DGP_estimate = estimateParams_RZ('RZDAT.xlsx','long','TVAR');

end


%% MONTE CARLO SIMULATIONS for CARs
%----------------------------------------------------------------
% Generate Data From DGP and compute CARs 
% 
%----------------------------------------------------------------
CAR_g = zeros(settings.est.IRF_hor, settings.est.n_state, length(delta_vals)); % 21x2x6
CAR_y = zeros(settings.est.IRF_hor, settings.est.n_state, length(delta_vals));
H = settings.est.IRF_hor; % 21
M = settings.simul.n_MC; 

rng(settings.simul.seed, 'twister');
for m = 1: M
    % Generate new dataset for each Monte Carlo iteration
    [synthetic_w, synthetic_S] = generateData_RZ(T_L,DGP_estimate,lag_type,0,false);
    % Temporary storage for impulse responses for this dataset for each regime i
    CAR_g_i = zeros(H, 2, length(delta_vals));
    CAR_y_i = zeros(H, 2, length(delta_vals));
    % Data Split here and loop over regimes
    for i_state =0:1
        idx = (synthetic_S == i_state);
        state_w = synthetic_w(idx,:);
        R_i=size(state_w,1);
 
        % Loop over shock sizes
        for d = 1:length(delta_vals)
            delta = delta_vals(d);
            num_skipped = 0; % Reset skipped count for each delta
            
            % Temporary storage for responses across histories for the current delta
            CAR_g_temp = zeros(H, R_i);
            CAR_y_temp = zeros(H, R_i);
            valid_index = 0; % Track valid cases separately

            for r = 1:R_i % Iterate over histories in regime i
                if r > 3
                   % Use the previous four lag information for histories to generate potential output, 4x3
                   Y_hist4 = state_w(r-3:r , :);
                else
                    % If t ≤ 3, throw away this data to compute CARs
                    num_skipped = num_skipped + 1;
                    continue
                end
                
                % Increment valid cases counter
                valid_index = valid_index + 1;

                % Path 1: With structural shock Y at r
                Y_with_shock = simulatePath(Y_hist4, DGP_estimate, delta, i_state);

                % Path 2: Without structural shock Y at t
                Y_without_shock = simulatePath(Y_hist4, DGP_estimate, 0, i_state);

                % Compute the difference
                CAR_g_temp(:, valid_index) = Y_with_shock(:, 2) - Y_without_shock(:, 2); % Government spending response
                CAR_y_temp(:, valid_index) = Y_with_shock(:, 3) - Y_without_shock(:, 3); % GDP response
            end
            
            % Average across valid histories for regime i
            valid_cases = R_i - num_skipped;
            if valid_cases > 0
                CAR_g_i(:, i_state+1, d) = sum(CAR_g_temp(:, 1:valid_cases), 2) / valid_cases;
                CAR_y_i(:, i_state+1, d) = sum(CAR_y_temp(:, 1:valid_cases), 2) / valid_cases;
            end
        end 
    end

    % Accumulate impulse responses across datasets (for each MC)
    CAR_g = CAR_g + CAR_g_i;
    CAR_y = CAR_y + CAR_y_i;
end

% Average across all datasets
CAR_g = CAR_g / M;
CAR_y = CAR_y / M;

% clear CAR_g_i CAR_g_temp CAR_y_i CAR_y_temp d delta num_skipped synthetic_S ... 
% synthetic_w t valid_cases valid_index Y_hist4 Y_without_shock Y_with_shock

