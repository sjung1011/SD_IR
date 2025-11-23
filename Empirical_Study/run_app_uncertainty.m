%% EMPIRICAL APPLICATION: MAIN FILE
% July 2025 ver
% Uncertainty shock
% Three transition processes, inflation-related
% State-dependent estimation: IRF, error bands with draws
% Include GIRF, substitution for SVAR
% VAR(12), LP(2)
%% HOUSEKEEPING
clc
clear all
close all

rng(1,'twister');
tic;
%% SET EXPERIMENT 
exper = 'Bloom'; % 
vol_type = 'actual'; % 'indicator', 'actual': go to mode_type 
estimand_type = 'IV'; % 'IV'(Vol) or 'Recursive'(FF) 
lag_type = 12; % for VARs
lag_type_lp =2;
mode_type = 5; % folder name:
               % 4 ('indicator'), 5 ('actual')
state_type = 'S3_75thr'; % 'S1_60_id','S1_30_id','S1_60_pc','S1_30_pc'
                         % 'S2_78_hh','S2_81_hh','S2_81_pro'
                         % 'S3_abs','S3_rel','S3_75thr' (10)
%% Settings
run(fullfile("DGP/Settings/check_mode.m"));
run(fullfile("Empirical_Study/shared_empirical.m"));

% Storage folder for results
save_pre = 'Results_Empirical'; % destination to store the results
save_suff = num2str(lag_type);
save_folder = fullfile(save_pre, save_mode_dir, strcat('lag', save_suff));

%% Load data
 % Note that 
        % Monthly data: July 1962-June 2008 (552 months) except for S2 
            % Michigan HH: Jan 1978 - June 2008 (366 months)
            % SPF: July 1981 - June 2008 (43:366 months)
load("Empirical_Study\Data_EMP.mat")
% Y: w_t =[proxy, \bar{w}], all observed data incuding proxy
if strcmp(vol_type,'indicator')
    Y_raw = table2array(Data.W(:,2:end));
else % actual
    VOLATBL = Data.W.VOLATBL;
    Y_mat = table2array(Data.W(:,3:end));
    Y_raw = [VOLATBL Y_mat];
    clear VOLATBL Y_mat 
end
% S
if strcmp(state_type,'S1_60_id')
    S_raw = Data.S1.S1_60_id;
elseif strcmp(state_type,'S1_30_id')
    S_raw = Data.S1.S1_30_id;
elseif strcmp(state_type,'S1_60_pc')
    S_raw = Data.S1.S1_60_pc;
elseif strcmp(state_type,'S1_30_pc')
    S_raw = Data.S1.S1_30_pc;
elseif strcmp(state_type,'S2_78_hh')
    S_raw = Data.S2.S_anchor78;
    Y_raw = Y_raw(187:end,:);
elseif strcmp(state_type,'S2_81_hh')
    S_raw = Data.S2.S_anchor81(43:end);
    Y_raw = Y_raw(229:end,:);
elseif strcmp(state_type,'S2_81_pro')
    S_raw = Data.S2.S_anchorPRO(43:end);
    Y_raw = Y_raw(229:end,:);
elseif strcmp(state_type,'S3_abs')
    S_raw = Data.S3.S_abs;
elseif strcmp(state_type,'S3_rel')
    S_raw = Data.S3.S_rel;
elseif strcmp(state_type,'S3_75thr')
    S_raw = Data.S3.S_75thr;
end
temp_hist_state0 = sum(S_raw(:) == 0);
temp_hist_state1 = sum(S_raw(:) == 1);
temp_hist_linear = size(Y_raw,1);
%% Placeholders
irfdraws_state0_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor,settings.est.n_btstrp); % 4 x 37 x 5e3
irfdraws_state1_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor,settings.est.n_btstrp); % 4 x 37 x 5e3
irfdraws_state0_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor,settings.est.n_btstrp); % 4 x 37 x 5e3
irfdraws_state1_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor,settings.est.n_btstrp); % 4 x 37 x 5e3
irfdraws_state0_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor,settings.est.n_btstrp); % 4 x 37 x 5e3
irfdraws_state1_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor,settings.est.n_btstrp); % 4 x 37 x 5e3
irfdraws_state0_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor,settings.est.n_btstrp); % 4 x 37 x 5e3
irfdraws_state1_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor,settings.est.n_btstrp); % 4 x 37 x 5e3
irf_state0_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); % 4 x H
irf_state1_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
irf_state0_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
irf_state1_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
irf_state0_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
irf_state1_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
irf_state0_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
irf_state1_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
gic_state0_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
gic_state1_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
gic_state0_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
gic_state1_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
gic_state0_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
gic_state1_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
gic_state0_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
gic_state1_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ssr_state0_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ssr_state1_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ssr_state0_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ssr_state1_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ssr_state0_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ssr_state1_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ssr_state0_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ssr_state1_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_lo_state0_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_lo_state1_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_lo_state0_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_lo_state1_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_lo_state0_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_lo_state1_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_lo_state0_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_lo_state1_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_hi_state0_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_hi_state1_vol_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_hi_state0_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_hi_state1_ff_empl =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_hi_state0_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_hi_state1_vol_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_hi_state0_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
ci_hi_state1_ff_ip =NaN(settings.est.nmethods_name,settings.est.IRF_hor); 
% directly receive tss and cCVA_w
% linear weight holders 
gic_linear_vol_empl =NaN(settings.est.nmethods_name_linear,settings.est.IRF_hor); 
gic_linear_ff_empl =NaN(settings.est.nmethods_name_linear,settings.est.IRF_hor); 
gic_linear_vol_ip =NaN(settings.est.nmethods_name_linear,settings.est.IRF_hor);  
gic_linear_ff_ip =NaN(settings.est.nmethods_name_linear,settings.est.IRF_hor); 
ssr_linear_vol_empl =NaN(settings.est.nmethods_name_linear,settings.est.IRF_hor); 
ssr_linear_ff_empl =NaN(settings.est.nmethods_name_linear,settings.est.IRF_hor); 
ssr_linear_vol_ip =NaN(settings.est.nmethods_name_linear,settings.est.IRF_hor); 
ssr_linear_ff_ip =NaN(settings.est.nmethods_name_linear,settings.est.IRF_hor); 
%
ei_svar = [];
ei_bvar = [];
ei_lp =[];
ei_lp_contmp =[];

%% EMPIRICAL ANALYSIS
%----------------------------------------------------------------
% Estimation for each vol type: Regime-sepecific weight
%----------------------------------------------------------------   
comb = {'state0_vol_empl', 'state1_vol_empl', ...
    'state0_ff_empl',  'state1_ff_empl',  ...
    'state0_vol_ip',   'state1_vol_ip',  ...
    'state0_ff_ip',    'state1_ff_ip' };

    for i_method = 1:length(settings.est.methods_name)
        switch settings.est.methods_name{i_method} % {'svar','bvar','lp','lp_contmp','girf'};
            case 'svar'
                [IRF_draws,IRF_hat,GIC,ei_svar,SSR,TSS] = SVAR_est_stateDep_emp(Y_raw, S_raw, settings, 0);
                % each struct has 8 variables.
                % IRF_draws.(): H × B; IRF_hat.(): H x 1; GIC, SSR, RSS.(): H x 1
                for c = 1:numel(comb) 
                    f = comb{c};                        % e.g. 'state0_vol_empl'
                    % --- IRF draws (H × B) -------------------------------------------
                    eval(['irfdraws_' f '(i_method,:,:) = IRF_draws.' f ';']);

                    % --- point IRF (H × 1) -------------------------------------------
                    eval(['irf_'       f '(i_method,:)  = IRF_hat.' f ';']);
                    % --- GIC, SSR,  (H × 1) ------------------------------------------
                    eval(['gic_'       f '(i_method,:)  = GIC.' f ';']);
                    eval(['ssr_'       f '(i_method,:)  = SSR.' f ';']);
                end
                tss_state0_vars_empl = TSS.state0_ff_empl; % H x 1
                tss_state1_vars_empl = TSS.state1_ff_empl;
                tss_state0_vars_ip = TSS.state0_ff_ip;
                tss_state1_vars_ip = TSS.state1_ff_ip;

                % boot_CI struct
                pct_lo = 16;          % choose the percentiles you want
                pct_hi = 84;

                ci_lo = struct();   % will hold H×1 lower bands
                ci_hi = struct();   % will hold H×1 upper bands

                fns = fieldnames(IRF_draws);   % e.g. 'state0_vol_emp', ...

                for k = 1:numel(fns)
                    fn        = fns{k};
                    draws     = IRF_draws.(fn);        % H × B
                    ci_lo.(fn) = prctile(draws, pct_lo, 2);   % H × 1
                    ci_hi.(fn) = prctile(draws, pct_hi, 2);   % H × 1
                end
                
                for c = 1:numel(comb)
                    f = comb{c};
                    eval(['ci_lo_' f '(i_method,:) = ci_lo.' f ';']);
                    eval(['ci_hi_' f '(i_method,:) = ci_hi.' f ';']);
                end
               
            case 'bvar'
                [IRF_draws,IRF_hat,GIC,ei_bvar,SSR, ~] ...
                    = BVAR_est_stateDep_emp(Y_raw, S_raw, settings);
                for c = 1:numel(comb) 
                    f = comb{c};                        % e.g. 'state0_vol_empl'
                    % --- IRF draws (H × B) -------------------------------------------
                    eval(['irfdraws_' f '(i_method,:,:) = IRF_draws.' f ';']);

                    % --- point IRF (H × 1) -------------------------------------------
                    eval(['irf_'       f '(i_method,:)  = IRF_hat.' f ';']);
                    % --- GIC, SSR,  (H × 1) ------------------------------------------
                    eval(['gic_'       f '(i_method,:)  = GIC.' f ';']);
                    eval(['ssr_'       f '(i_method,:)  = SSR.' f ';']);
                end
                % boot_CI struct
                pct_lo = 16;          % choose the percentiles you want
                pct_hi = 84;

                ci_lo = struct();   % will hold H×1 lower bands
                ci_hi = struct();   % will hold H×1 upper bands

                fns = fieldnames(IRF_draws);   % e.g. 'state0_vol_emp', ...

                for k = 1:numel(fns)
                    fn        = fns{k};
                    draws     = IRF_draws.(fn);        % H × B
                    ci_lo.(fn) = prctile(draws, pct_lo, 2);   % H × 1
                    ci_hi.(fn) = prctile(draws, pct_hi, 2);   % H × 1
                end
                
                for c = 1:numel(comb)
                    f = comb{c};
                    eval(['ci_lo_' f '(i_method,:) = ci_lo.' f ';']);
                    eval(['ci_hi_' f '(i_method,:) = ci_hi.' f ';']);
                end
                % 
            case 'lp'
                [IRF_draws,IRF_hat,GIC,ei_lp,SSR,TSS]... 
                =LP_est_stateDep_emp(Y_raw, S_raw, settings);
                for c = 1:numel(comb) 
                    f = comb{c};                        % e.g. 'state0_vol_empl'
                    % --- IRF draws (H × B) -------------------------------------------
                    eval(['irfdraws_' f '(i_method,:,:) = IRF_draws.' f ';']);

                    % --- point IRF (H × 1) -------------------------------------------
                    eval(['irf_'       f '(i_method,:)  = IRF_hat.' f ';']);
                    % --- GIC, SSR,  (H × 1) ------------------------------------------
                    eval(['gic_'       f '(i_method,:)  = GIC.' f ';']);
                    eval(['ssr_'       f '(i_method,:)  = SSR.' f ';']);
                end
                tss_state0_lps_empl = TSS.state0_ff_empl; % H x 1
                tss_state1_lps_empl = TSS.state1_ff_empl;
                tss_state0_lps_ip = TSS.state0_ff_ip;
                tss_state1_lps_ip = TSS.state1_ff_ip;
                % boot_CI struct
                pct_lo = 16;          % choose the percentiles you want
                pct_hi = 84;

                ci_lo = struct();   % will hold H×1 lower bands
                ci_hi = struct();   % will hold H×1 upper bands

                fns = fieldnames(IRF_draws);   % e.g. 'state0_vol_emp', ...

                for k = 1:numel(fns)
                    fn        = fns{k};
                    draws     = IRF_draws.(fn);        % H × B
                    ci_lo.(fn) = prctile(draws, pct_lo, 2);   % H × 1
                    ci_hi.(fn) = prctile(draws, pct_hi, 2);   % H × 1
                end
                
                for c = 1:numel(comb)
                    f = comb{c};
                    eval(['ci_lo_' f '(i_method,:) = ci_lo.' f ';']);
                    eval(['ci_hi_' f '(i_method,:) = ci_hi.' f ';']);
                end
                %
            case 'lp_contmp'
                [IRF_draws,IRF_hat,GIC,ei_lp_contmp,SSR,~]... 
                =LP_est_stateDep_emp_contmp(Y_raw, S_raw, settings);
                for c = 1:numel(comb) 
                    f = comb{c};                        % e.g. 'state0_vol_empl'
                    % --- IRF draws (H × B) -------------------------------------------
                    eval(['irfdraws_' f '(i_method,:,:) = IRF_draws.' f ';']);

                    % --- point IRF (H × 1) -------------------------------------------
                    eval(['irf_'       f '(i_method,:)  = IRF_hat.' f ';']);
                    % --- GIC, SSR,  (H × 1) ------------------------------------------
                    eval(['gic_'       f '(i_method,:)  = GIC.' f ';']);
                    eval(['ssr_'       f '(i_method,:)  = SSR.' f ';']);
                end
                % boot_CI struct
                pct_lo = 16;          % choose the percentiles you want
                pct_hi = 84;

                ci_lo = struct();   % will hold H×1 lower bands
                ci_hi = struct();   % will hold H×1 upper bands

                fns = fieldnames(IRF_draws);   % e.g. 'state0_vol_emp', ...

                for k = 1:numel(fns)
                    fn        = fns{k};
                    draws     = IRF_draws.(fn);        % H × B
                    ci_lo.(fn) = prctile(draws, pct_lo, 2);   % H × 1
                    ci_hi.(fn) = prctile(draws, pct_hi, 2);   % H × 1
                end
                
                for c = 1:numel(comb)
                    f = comb{c};
                    eval(['ci_lo_' f '(i_method,:) = ci_lo.' f ';']);
                    eval(['ci_hi_' f '(i_method,:) = ci_hi.' f ';']);
                end
             % substitution of svar
            case 'girf'
                [IRF_draws, IRF_hat] = GIRF_VAR_boot_stateDep(Y_raw, S_raw, settings);
                for c = 1:numel(comb) 
                    f = comb{c};                        % e.g. 'state0_vol_empl'
                    % --- IRF draws (H × B) -------------------------------------------
                    eval(['irfdraws_' f '(i_method,:,:) = IRF_draws.' f ';']);

                    % --- gIRF (H × 1) -------------------------------------------
                    eval(['irf_'       f '(i_method,:)  = IRF_hat.' f ';']);
                end
                % boot_CI struct
                pct_lo = 16;          % 68% ci
                pct_hi = 84;

                ci_lo = struct();   % will hold H×1 lower bands
                ci_hi = struct();   % will hold H×1 upper bands

                fns = fieldnames(IRF_draws);   % e.g. 'state0_vol_emp', ...

                for k = 1:numel(fns)
                    fn        = fns{k};
                    draws     = IRF_draws.(fn);        % H × B
                    ci_lo.(fn) = prctile(draws, pct_lo, 2);   % H × 1
                    ci_hi.(fn) = prctile(draws, pct_hi, 2);   % H × 1
                end
                
                for c = 1:numel(comb)
                    f = comb{c};
                    eval(['ci_lo_' f '(i_method,:) = ci_lo.' f ';']);
                    eval(['ci_hi_' f '(i_method,:) = ci_hi.' f ';']);
                end
        end
    end 
[temp_CVA_w_vars] = CVA_est_var_stateDep_emp(Y_raw, S_raw, settings,ei_svar,ei_bvar);% each filed H x 2
[temp_CVA_w_lps] = CVA_est_lp_stateDep_emp(Y_raw, S_raw, settings,ei_lp,ei_lp_contmp);

%% --------------------------------------------------------------
% Estimation for Regime-constant weight
%----------------------------------------------------------------
comb_linear = {'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip' };

    for i_method = 1:length(settings.est.methods_name_linear)
        switch settings.est.methods_name_linear{i_method} % {'svar','bvar','lp','lp_contmp'};
            case 'svar'
                [~,~,GIC,ei_svar,SSR,TSS] = SVAR_est_linear_emp(Y_raw, settings, 0);
                % IRF_draws.(): H × B; IRF_hat.(): H x 1; GIC, SSR, RSS.(): H x 1
                for c = 1:numel(comb_linear) 
                    f = comb_linear{c};                        % e.g. 'state0_vol_empl'
                  
                    % --- GIC, SSR,  (H × 1) ------------------------------------------
                    eval(['gic_'       f '(i_method,:)  = GIC.' f ';']);
                    eval(['ssr_'       f '(i_method,:)  = SSR.' f ';']);
                end
                tss_linear_vars_empl = TSS.linear_ff_empl; % H x 1
                tss_linear_vars_ip = TSS.linear_ff_ip;
                
            case 'bvar'
                [~,~,GIC,ei_bvar,SSR, ~] ...
                    = BVAR_est_linear_emp(Y_raw,settings);
                for c = 1:numel(comb_linear) 
                    f = comb_linear{c};                       % e.g. 'state0_vol_empl'
                    
                    % --- GIC, SSR,  (H × 1) ------------------------------------------
                    eval(['gic_'       f '(i_method,:)  = GIC.' f ';']);
                    eval(['ssr_'       f '(i_method,:)  = SSR.' f ';']);
                end
                % 
            case 'lp'
                [~,~,GIC,ei_lp,SSR,TSS]... 
                =LP_est_linear_emp(Y_raw,settings);
                for c = 1:numel(comb_linear) 
                    f = comb_linear{c};                          % e.g. 'state0_vol_empl'
                    
                    % --- GIC, SSR,  (H × 1) ------------------------------------------
                    eval(['gic_'       f '(i_method,:)  = GIC.' f ';']);
                    eval(['ssr_'       f '(i_method,:)  = SSR.' f ';']);
                end
                tss_linear_lps_empl = TSS.linear_ff_empl; % H x 1
                tss_linear_lps_ip = TSS.linear_ff_ip;
                
                %
            case 'lp_contmp'
                [~,~,GIC,ei_lp_contmp,SSR,~]... 
                =LP_est_linear_emp_contmp(Y_raw, settings);
                for c = 1:numel(comb_linear) 
                    f = comb_linear{c};                          % e.g. 'state0_vol_empl'
                    
                    % --- GIC, SSR,  (H × 1) ------------------------------------------
                    eval(['gic_'       f '(i_method,:)  = GIC.' f ';']);
                    eval(['ssr_'       f '(i_method,:)  = SSR.' f ';']);
                end
        end
    end 
[temp_CVA_w_vars_linear] = CVA_est_var_linear_emp(Y_raw,settings,ei_svar,ei_bvar);% each filed H x 2
[temp_CVA_w_lps_linear] = CVA_est_lp_linear_emp(Y_raw, settings,ei_lp,ei_lp_contmp);



%% SUMMARIZE RESULTS
for i_method = 1:settings.est.nmethods_name
    thisMethod = settings.est.methods_name{i_method};
    for c = 1:numel(comb)
        f = comb{c};
        %
        tmp = eval(['irfdraws_' f ...
               '(i_method, settings.est.IRF_select, :)']);
        irfpaths.(thisMethod).(state_type).(f) =  permute(tmp,[2 3 1]);
        %
        tmp = eval(['irf_' f ...
               '(i_method, settings.est.IRF_select, :)']);
        irf.(thisMethod).(state_type).(f) = permute(tmp,[2 1]);
        %
        tmp = eval(['gic_' f ...
               '(i_method, settings.est.IRF_select, :)']);
        gic.(thisMethod).(state_type).(f) = permute(tmp,[2 1]);
        %
        tmp = eval(['ssr_' f ...
               '(i_method, settings.est.IRF_select, :)']);
        ssr.(thisMethod).(state_type).(f) = permute(tmp,[2 1]);
        %
        tmp = eval(['ci_lo_' f ...
               '(i_method, settings.est.IRF_select, :)']);
        band_lo.(thisMethod).(state_type).(f) = permute(tmp,[2 1]);
        %
        tmp = eval(['ci_hi_' f ...
               '(i_method, settings.est.IRF_select, :)']);
        band_hi.(thisMethod).(state_type).(f) = permute(tmp,[2 1]);
    end
end
%
comb = {'state0_lps_empl', 'state0_vars_empl', ...
    'state1_lps_empl',  'state1_vars_empl',  ...
    'state0_lps_ip',   'state0_vars_ip',  ...
    'state1_lps_ip',    'state1_vars_ip' };
for c = 1:numel(comb)
    f = comb{c};
    tmp = eval(['tss_' f ]);
    tss.(state_type).(f) = tmp;
end
CVA_w_vars.(state_type) = temp_CVA_w_vars;
CVA_w_lps.(state_type) = temp_CVA_w_lps;

% linear weight -------------------------------------------------
for i_method = 1:settings.est.nmethods_name_linear
    thisMethod = settings.est.methods_name_linear{i_method};
    for c = 1:numel(comb_linear)
        f = comb_linear{c};
        %
        tmp = eval(['gic_' f ...
               '(i_method, settings.est.IRF_select, :)']);
        gic.(thisMethod).(state_type).(f) = permute(tmp,[2 1]);
        %
        tmp = eval(['ssr_' f ...
               '(i_method, settings.est.IRF_select, :)']);
        ssr.(thisMethod).(state_type).(f) = permute(tmp,[2 1]);
     end
end
for c = 1:numel(comb_linear)
    f = comb_linear{c};                       
    tmp = temp_CVA_w_vars_linear.(f);           
    CVA_w_vars.(state_type).(f) = tmp;
    tmp = temp_CVA_w_lps_linear.(f);
    CVA_w_lps.(state_type).(f)  = tmp;
end
%
comb_linear = {'linear_lps_empl', 'linear_vars_empl', ...
        'linear_lps_ip',   'linear_vars_ip'};

for c = 1:numel(comb_linear)
    f = comb_linear{c};
    tmp = eval(['tss_' f ]);
    tss.(state_type).(f) = tmp;
end
% 
hist.(state_type).state0 = temp_hist_state0;
hist.(state_type).state1 = temp_hist_state1;
hist.(state_type).linear = temp_hist_linear;
%%
clear tss_* ssr_* gic_* irfdraws_* ci_lo_* ci_hi_* temp_* ei_* irf_*
clear i_method thisMethod tmp c comb f fn fns k comb_linear
clear Data S1Descr S2Descr S3Descr S_raw WDescr Y_raw
clear ci_lo ci_hi IRF_draws IRF_hat GIC SSR TSS draws
%% Export Results

mkdir(save_folder);
save(fullfile(save_folder, strcat(exper,'_', state_type, '_', vol_type)), ...
    'settings','lag_type', 'pct_hi', 'pct_lo', 'lag_type_lp',...
    'vol_type','state_type',...
    'irfpaths','irf','gic','ssr','band_lo','band_hi','tss', 'hist', ...
    'CVA_w_vars', 'CVA_w_lps',...
    '-v7.3'); % save results

clear save_folder save_pre save_mode_dir save_suff mode_list mode_type
toc;