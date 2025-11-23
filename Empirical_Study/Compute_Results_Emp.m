%% STATE-DEPENDENT EMPIRICAL STUDY: Compute and Save linear and SD Results 
% SEOJIN JUNG
% This version: Aug.2025
% IRFs per regime, hs, respV
    % 5e3 MCs results
    % Regime specific weight  
    % No CVA_w_lps (sd): S2_81_hh_indicator
%%
clc
clear all
close all
%% Collect data and Create MAT-file
filename = 'Combine_Results_Emp.mat';
save(filename);
clear filename;
%% 
run(fullfile("Empirical_Study/shared_empirical.m"));
%% -------------------------------------------
% measuring the shock size
% --------------------------------------------
load("Empirical_Study\Data_EMP.mat","Data")

% Compute standard deviations
sigma.monetary = std(Data.W.cffr, 'omitnan');
sigma.uncert_ind   = std(Data.W.clean, 'omitnan');
sigma.uncert   = std(Data.W.VOLATBL, 'omitnan');

% Normalize shocks by their std = \delta_t/sigma_1
stdUnits.monetary = Data.W.cffr ./sigma.monetary;
stdUnits.uncert_ind   = Data.W.clean ./ sigma.uncert_ind;
stdUnits.uncert = Data.W.VOLATBL ./ sigma.uncert;

% average unit shock of them 
stdUnits.mean_monetary = mean(stdUnits.monetary,1); % 0.0012
stdUnits.mean_uncert_ind = mean(stdUnits.uncert_ind,1); % 0.1781
stdUnits.mean_uncert = mean(stdUnits.uncert,1); % 3.0228

% T_62, T_78, T_81
T_62 = size(Data.S1,1);
T_78 = size(Data.S2,1);
T_81 = T_78-43+1;
%%
save("Combine_Results_Emp.mat",'sigma','stdUnits','T_62','T_78','T_81')

%% ------------------------------------------------------
% 2nd shock = actual
% -------------------------------------------------------
% 1. Good and bad inflation regimes 
% 1-1. the rolling window (60 months) correlation between inflation (yoy) and
% consumption growth (index)
% S1_60_id_actual
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\actual\lag12\Bloom_S1_60_id_actual.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S1_60_id_actual.(names{i}).irf.(f) = irfpaths.(names{i}).S1_60_id.(f);
        S1_60_id_actual.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S1_60_id.(f),2);
        S1_60_id_actual.(names{i}).band_hi.(f) = band_hi.(names{i}).S1_60_id.(f);
        S1_60_id_actual.(names{i}).band_lo.(f) = band_lo.(names{i}).S1_60_id.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S1_60_id.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S1_60_id_actual.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_60_id.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S1_60_id_actual.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_60_id.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S1_60_id_actual.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S1_60_id.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S1_60_id_actual.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S1_60_id.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S1_60_id_actual.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S1_60_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S1_60_id_actual.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_actual.GMA_lps.irf.(f) = wirf_val;% H x B
    S1_60_id_actual.GMA_lps.mu.(f) = mean(S1_60_id_actual.GMA_lps.irf.(f), 2, "omitmissing");
    S1_60_id_actual.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_actual.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S1_60_id_actual.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_actual.CVA_lps.irf.(f) = wirf_val;% H x B
    S1_60_id_actual.CVA_lps.mu.(f) = mean(S1_60_id_actual.CVA_lps.irf.(f), 2, "omitmissing");
    S1_60_id_actual.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_actual.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S1_60_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_60_id_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_actual.GMA_vars.irf.(f) = wirf_val;% H x B
    S1_60_id_actual.GMA_vars.mu.(f) = mean(S1_60_id_actual.GMA_vars.irf.(f), 2, "omitmissing");
    S1_60_id_actual.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_actual.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_60_id_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_actual.CVA_vars.irf.(f) = wirf_val;
    S1_60_id_actual.CVA_vars.mu.(f) = mean(S1_60_id_actual.CVA_vars.irf.(f), 2, "omitmissing");
    S1_60_id_actual.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_actual.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S1_60_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_60_id_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_actual.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S1_60_id_actual.GMA_vars_girf.mu.(f) = mean(S1_60_id_actual.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S1_60_id_actual.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_actual.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_60_id_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_actual.CVA_vars_girf.irf.(f) = wirf_val;
    S1_60_id_actual.CVA_vars_girf.mu.(f) = mean(S1_60_id_actual.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S1_60_id_actual.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_actual.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S1_60_id.state0_vol_empl = tss.S1_60_id.state0_lps_empl;
    tss.S1_60_id.state0_ff_empl = tss.S1_60_id.state0_lps_empl;
    tss.S1_60_id.state1_vol_empl = tss.S1_60_id.state1_lps_empl;
    tss.S1_60_id.state1_ff_empl = tss.S1_60_id.state1_lps_empl;
    tss.S1_60_id.state0_vol_ip = tss.S1_60_id.state0_lps_ip;
    tss.S1_60_id.state0_ff_ip = tss.S1_60_id.state0_lps_ip;
    tss.S1_60_id.state1_vol_ip = tss.S1_60_id.state1_lps_ip;
    tss.S1_60_id.state1_ff_ip = tss.S1_60_id.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S1_60_id.(f) ./ tss.S1_60_id.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S1_60_id.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S1_60_id.state0_vol_empl = tss.S1_60_id.state0_vars_empl;
    tss.S1_60_id.state0_ff_empl = tss.S1_60_id.state0_vars_empl;
    tss.S1_60_id.state1_vol_empl = tss.S1_60_id.state1_vars_empl;
    tss.S1_60_id.state1_ff_empl = tss.S1_60_id.state1_vars_empl;
    tss.S1_60_id.state0_vol_ip = tss.S1_60_id.state0_vars_ip;
    tss.S1_60_id.state0_ff_ip = tss.S1_60_id.state0_vars_ip;
    tss.S1_60_id.state1_vol_ip = tss.S1_60_id.state1_vars_ip;
    tss.S1_60_id.state1_ff_ip = tss.S1_60_id.state1_vars_ip;
    tss.S1_60_id.linear_vol_empl = tss.S1_60_id.linear_vars_empl;
    tss.S1_60_id.linear_vol_ip = tss.S1_60_id.linear_vars_ip;
    tss.S1_60_id.linear_ff_empl = tss.S1_60_id.linear_vars_empl;
    tss.S1_60_id.linear_ff_ip = tss.S1_60_id.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S1_60_id.(f) ./ tss.S1_60_id.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S1_60_id.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S1_60_id_actual.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S1_60_id_actual.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S1_60_id_actual.aLP_RS.alpha.(f) .* S1_60_id_actual.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S1_60_id_actual.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_60_id_actual.GMA_vars.weight.(f)'; % H x 2
    S1_60_id_actual.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S1_60_id_actual.aLP_RS.alpha.(f) .* S1_60_id_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_60_id_actual.CVA_vars.weight.(f)';
    S1_60_id_actual.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S1_60_id_actual.aLP_MSPE.alpha.(f) .* S1_60_id_actual.GMA_lps.weight.(f)';
    alpha_var = 1 - S1_60_id_actual.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_60_id_actual.GMA_vars.weight.(f)';
    S1_60_id_actual.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S1_60_id_actual.aLP_MSPE.alpha.(f) .* S1_60_id_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_60_id_actual.CVA_vars.weight.(f)';
    S1_60_id_actual.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S1_60_id_actual.aLP_RS_GMA.irf.(f) = S1_60_id_actual.aLP_RS.alpha.(f) .* S1_60_id_actual.GMA_lps.irf.(f)...
    + (1-S1_60_id_actual.aLP_RS.alpha.(f)) .* S1_60_id_actual.GMA_vars.irf.(f); 
    S1_60_id_actual.aLP_RS_GMA.mu.(f) = mean(S1_60_id_actual.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S1_60_id_actual.aLP_RS_GMA.band_lo.(f) = prctile(S1_60_id_actual.aLP_RS_GMA.irf.(f),pct_lo,2);
    S1_60_id_actual.aLP_RS_GMA.band_hi.(f) = prctile(S1_60_id_actual.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S1_60_id_actual.aLP_RS_CVA.irf.(f) = S1_60_id_actual.aLP_RS.alpha.(f) .* S1_60_id_actual.CVA_lps.irf.(f)...
    + (1-S1_60_id_actual.aLP_RS.alpha.(f)) .* S1_60_id_actual.CVA_vars.irf.(f); 
    S1_60_id_actual.aLP_RS_CVA.mu.(f) = mean(S1_60_id_actual.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S1_60_id_actual.aLP_RS_CVA.band_lo.(f) = prctile(S1_60_id_actual.aLP_RS_CVA.irf.(f),pct_lo,2);
    S1_60_id_actual.aLP_RS_CVA.band_hi.(f) = prctile(S1_60_id_actual.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S1_60_id_actual.aLP_MSPE_GMA.irf.(f) = S1_60_id_actual.aLP_MSPE.alpha.(f) .* S1_60_id_actual.GMA_lps.irf.(f)...
    + (1-S1_60_id_actual.aLP_MSPE.alpha.(f)) .* S1_60_id_actual.GMA_vars.irf.(f); 
    S1_60_id_actual.aLP_MSPE_GMA.mu.(f) = mean(S1_60_id_actual.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S1_60_id_actual.aLP_MSPE_GMA.band_lo.(f) = prctile(S1_60_id_actual.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S1_60_id_actual.aLP_MSPE_GMA.band_hi.(f) = prctile(S1_60_id_actual.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S1_60_id_actual.aLP_MSPE_CVA.irf.(f) = S1_60_id_actual.aLP_MSPE.alpha.(f) .* S1_60_id_actual.CVA_lps.irf.(f)...
    + (1-S1_60_id_actual.aLP_MSPE.alpha.(f)) .* S1_60_id_actual.CVA_vars.irf.(f); 
    S1_60_id_actual.aLP_MSPE_CVA.mu.(f) = mean(S1_60_id_actual.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S1_60_id_actual.aLP_MSPE_CVA.band_lo.(f) = prctile(S1_60_id_actual.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S1_60_id_actual.aLP_MSPE_CVA.band_hi.(f) = prctile(S1_60_id_actual.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S1_60_id_actual.aLP_RS_GMA_girf.irf.(f) = S1_60_id_actual.aLP_RS.alpha.(f) .* S1_60_id_actual.GMA_lps.irf.(f)...
    + (1-S1_60_id_actual.aLP_RS.alpha.(f)) .* S1_60_id_actual.GMA_vars_girf.irf.(f); 
    S1_60_id_actual.aLP_RS_GMA_girf.mu.(f) = mean(S1_60_id_actual.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S1_60_id_actual.aLP_RS_GMA_girf.band_lo.(f) = prctile(S1_60_id_actual.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S1_60_id_actual.aLP_RS_GMA_girf.band_hi.(f) = prctile(S1_60_id_actual.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_60_id_actual.aLP_RS_CVA_girf.irf.(f) = S1_60_id_actual.aLP_RS.alpha.(f) .* S1_60_id_actual.CVA_lps.irf.(f)...
    + (1-S1_60_id_actual.aLP_RS.alpha.(f)) .* S1_60_id_actual.CVA_vars_girf.irf.(f); 
    S1_60_id_actual.aLP_RS_CVA_girf.mu.(f) = mean(S1_60_id_actual.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S1_60_id_actual.aLP_RS_CVA_girf.band_lo.(f) = prctile(S1_60_id_actual.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S1_60_id_actual.aLP_RS_CVA_girf.band_hi.(f) = prctile(S1_60_id_actual.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S1_60_id_actual.aLP_MSPE_GMA_girf.irf.(f) = S1_60_id_actual.aLP_MSPE.alpha.(f) .* S1_60_id_actual.GMA_lps.irf.(f)...
    + (1-S1_60_id_actual.aLP_MSPE.alpha.(f)) .* S1_60_id_actual.GMA_vars_girf.irf.(f); 
    S1_60_id_actual.aLP_MSPE_GMA_girf.mu.(f) = mean(S1_60_id_actual.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S1_60_id_actual.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S1_60_id_actual.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S1_60_id_actual.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S1_60_id_actual.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_60_id_actual.aLP_MSPE_CVA_girf.irf.(f) = S1_60_id_actual.aLP_MSPE.alpha.(f) .* S1_60_id_actual.CVA_lps.irf.(f)...
    + (1-S1_60_id_actual.aLP_MSPE.alpha.(f)) .* S1_60_id_actual.CVA_vars_girf.irf.(f); 
    S1_60_id_actual.aLP_MSPE_CVA_girf.mu.(f) = mean(S1_60_id_actual.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S1_60_id_actual.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S1_60_id_actual.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S1_60_id_actual.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S1_60_id_actual.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S1_60_id_actual','-append')
%% ------------------------------------------------------- 
% 1-2. the rolling window (60 months) correlation between inflation (yoy) and
% consumption growth (per capita)
% S1_60_pc_actual
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\actual\lag12\Bloom_S1_60_pc_actual.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S1_60_pc_actual.(names{i}).irf.(f) = irfpaths.(names{i}).S1_60_pc.(f);
        S1_60_pc_actual.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S1_60_pc.(f),2);
        S1_60_pc_actual.(names{i}).band_hi.(f) = band_hi.(names{i}).S1_60_pc.(f);
        S1_60_pc_actual.(names{i}).band_lo.(f) = band_lo.(names{i}).S1_60_pc.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S1_60_pc.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S1_60_pc_actual.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_60_pc.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S1_60_pc_actual.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_60_pc.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S1_60_pc_actual.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S1_60_pc.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S1_60_pc_actual.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S1_60_pc.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S1_60_pc_actual.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S1_60_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S1_60_pc_actual.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_actual.GMA_lps.irf.(f) = wirf_val;% H x B
    S1_60_pc_actual.GMA_lps.mu.(f) = mean(S1_60_pc_actual.GMA_lps.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_actual.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S1_60_pc_actual.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_actual.CVA_lps.irf.(f) = wirf_val;% H x B
    S1_60_pc_actual.CVA_lps.mu.(f) = mean(S1_60_pc_actual.CVA_lps.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_actual.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S1_60_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_60_pc_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_actual.GMA_vars.irf.(f) = wirf_val;% H x B
    S1_60_pc_actual.GMA_vars.mu.(f) = mean(S1_60_pc_actual.GMA_vars.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_actual.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_60_pc_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_actual.CVA_vars.irf.(f) = wirf_val;
    S1_60_pc_actual.CVA_vars.mu.(f) = mean(S1_60_pc_actual.CVA_vars.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_actual.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S1_60_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_60_pc_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_actual.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S1_60_pc_actual.GMA_vars_girf.mu.(f) = mean(S1_60_pc_actual.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_actual.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_60_pc_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_actual.CVA_vars_girf.irf.(f) = wirf_val;
    S1_60_pc_actual.CVA_vars_girf.mu.(f) = mean(S1_60_pc_actual.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_actual.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S1_60_pc.state0_vol_empl = tss.S1_60_pc.state0_lps_empl;
    tss.S1_60_pc.state0_ff_empl = tss.S1_60_pc.state0_lps_empl;
    tss.S1_60_pc.state1_vol_empl = tss.S1_60_pc.state1_lps_empl;
    tss.S1_60_pc.state1_ff_empl = tss.S1_60_pc.state1_lps_empl;
    tss.S1_60_pc.state0_vol_ip = tss.S1_60_pc.state0_lps_ip;
    tss.S1_60_pc.state0_ff_ip = tss.S1_60_pc.state0_lps_ip;
    tss.S1_60_pc.state1_vol_ip = tss.S1_60_pc.state1_lps_ip;
    tss.S1_60_pc.state1_ff_ip = tss.S1_60_pc.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S1_60_pc.(f) ./ tss.S1_60_pc.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S1_60_pc.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S1_60_pc.state0_vol_empl = tss.S1_60_pc.state0_vars_empl;
    tss.S1_60_pc.state0_ff_empl = tss.S1_60_pc.state0_vars_empl;
    tss.S1_60_pc.state1_vol_empl = tss.S1_60_pc.state1_vars_empl;
    tss.S1_60_pc.state1_ff_empl = tss.S1_60_pc.state1_vars_empl;
    tss.S1_60_pc.state0_vol_ip = tss.S1_60_pc.state0_vars_ip;
    tss.S1_60_pc.state0_ff_ip = tss.S1_60_pc.state0_vars_ip;
    tss.S1_60_pc.state1_vol_ip = tss.S1_60_pc.state1_vars_ip;
    tss.S1_60_pc.state1_ff_ip = tss.S1_60_pc.state1_vars_ip;
    tss.S1_60_pc.linear_vol_empl = tss.S1_60_pc.linear_vars_empl;
    tss.S1_60_pc.linear_vol_ip = tss.S1_60_pc.linear_vars_ip;
    tss.S1_60_pc.linear_ff_empl = tss.S1_60_pc.linear_vars_empl;
    tss.S1_60_pc.linear_ff_ip = tss.S1_60_pc.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S1_60_pc.(f) ./ tss.S1_60_pc.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S1_60_pc.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S1_60_pc_actual.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S1_60_pc_actual.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S1_60_pc_actual.aLP_RS.alpha.(f) .* S1_60_pc_actual.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S1_60_pc_actual.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_60_pc_actual.GMA_vars.weight.(f)'; % H x 2
    S1_60_pc_actual.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S1_60_pc_actual.aLP_RS.alpha.(f) .* S1_60_pc_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_60_pc_actual.CVA_vars.weight.(f)';
    S1_60_pc_actual.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S1_60_pc_actual.aLP_MSPE.alpha.(f) .* S1_60_pc_actual.GMA_lps.weight.(f)';
    alpha_var = 1 - S1_60_pc_actual.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_60_pc_actual.GMA_vars.weight.(f)';
    S1_60_pc_actual.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S1_60_pc_actual.aLP_MSPE.alpha.(f) .* S1_60_pc_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_60_pc_actual.CVA_vars.weight.(f)';
    S1_60_pc_actual.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S1_60_pc_actual.aLP_RS_GMA.irf.(f) = S1_60_pc_actual.aLP_RS.alpha.(f) .* S1_60_pc_actual.GMA_lps.irf.(f)...
    + (1-S1_60_pc_actual.aLP_RS.alpha.(f)) .* S1_60_pc_actual.GMA_vars.irf.(f); 
    S1_60_pc_actual.aLP_RS_GMA.mu.(f) = mean(S1_60_pc_actual.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.aLP_RS_GMA.band_lo.(f) = prctile(S1_60_pc_actual.aLP_RS_GMA.irf.(f),pct_lo,2);
    S1_60_pc_actual.aLP_RS_GMA.band_hi.(f) = prctile(S1_60_pc_actual.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S1_60_pc_actual.aLP_RS_CVA.irf.(f) = S1_60_pc_actual.aLP_RS.alpha.(f) .* S1_60_pc_actual.CVA_lps.irf.(f)...
    + (1-S1_60_pc_actual.aLP_RS.alpha.(f)) .* S1_60_pc_actual.CVA_vars.irf.(f); 
    S1_60_pc_actual.aLP_RS_CVA.mu.(f) = mean(S1_60_pc_actual.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.aLP_RS_CVA.band_lo.(f) = prctile(S1_60_pc_actual.aLP_RS_CVA.irf.(f),pct_lo,2);
    S1_60_pc_actual.aLP_RS_CVA.band_hi.(f) = prctile(S1_60_pc_actual.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S1_60_pc_actual.aLP_MSPE_GMA.irf.(f) = S1_60_pc_actual.aLP_MSPE.alpha.(f) .* S1_60_pc_actual.GMA_lps.irf.(f)...
    + (1-S1_60_pc_actual.aLP_MSPE.alpha.(f)) .* S1_60_pc_actual.GMA_vars.irf.(f); 
    S1_60_pc_actual.aLP_MSPE_GMA.mu.(f) = mean(S1_60_pc_actual.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.aLP_MSPE_GMA.band_lo.(f) = prctile(S1_60_pc_actual.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S1_60_pc_actual.aLP_MSPE_GMA.band_hi.(f) = prctile(S1_60_pc_actual.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S1_60_pc_actual.aLP_MSPE_CVA.irf.(f) = S1_60_pc_actual.aLP_MSPE.alpha.(f) .* S1_60_pc_actual.CVA_lps.irf.(f)...
    + (1-S1_60_pc_actual.aLP_MSPE.alpha.(f)) .* S1_60_pc_actual.CVA_vars.irf.(f); 
    S1_60_pc_actual.aLP_MSPE_CVA.mu.(f) = mean(S1_60_pc_actual.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.aLP_MSPE_CVA.band_lo.(f) = prctile(S1_60_pc_actual.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S1_60_pc_actual.aLP_MSPE_CVA.band_hi.(f) = prctile(S1_60_pc_actual.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S1_60_pc_actual.aLP_RS_GMA_girf.irf.(f) = S1_60_pc_actual.aLP_RS.alpha.(f) .* S1_60_pc_actual.GMA_lps.irf.(f)...
    + (1-S1_60_pc_actual.aLP_RS.alpha.(f)) .* S1_60_pc_actual.GMA_vars_girf.irf.(f); 
    S1_60_pc_actual.aLP_RS_GMA_girf.mu.(f) = mean(S1_60_pc_actual.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.aLP_RS_GMA_girf.band_lo.(f) = prctile(S1_60_pc_actual.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S1_60_pc_actual.aLP_RS_GMA_girf.band_hi.(f) = prctile(S1_60_pc_actual.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_60_pc_actual.aLP_RS_CVA_girf.irf.(f) = S1_60_pc_actual.aLP_RS.alpha.(f) .* S1_60_pc_actual.CVA_lps.irf.(f)...
    + (1-S1_60_pc_actual.aLP_RS.alpha.(f)) .* S1_60_pc_actual.CVA_vars_girf.irf.(f); 
    S1_60_pc_actual.aLP_RS_CVA_girf.mu.(f) = mean(S1_60_pc_actual.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.aLP_RS_CVA_girf.band_lo.(f) = prctile(S1_60_pc_actual.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S1_60_pc_actual.aLP_RS_CVA_girf.band_hi.(f) = prctile(S1_60_pc_actual.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S1_60_pc_actual.aLP_MSPE_GMA_girf.irf.(f) = S1_60_pc_actual.aLP_MSPE.alpha.(f) .* S1_60_pc_actual.GMA_lps.irf.(f)...
    + (1-S1_60_pc_actual.aLP_MSPE.alpha.(f)) .* S1_60_pc_actual.GMA_vars_girf.irf.(f); 
    S1_60_pc_actual.aLP_MSPE_GMA_girf.mu.(f) = mean(S1_60_pc_actual.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S1_60_pc_actual.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S1_60_pc_actual.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S1_60_pc_actual.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_60_pc_actual.aLP_MSPE_CVA_girf.irf.(f) = S1_60_pc_actual.aLP_MSPE.alpha.(f) .* S1_60_pc_actual.CVA_lps.irf.(f)...
    + (1-S1_60_pc_actual.aLP_MSPE.alpha.(f)) .* S1_60_pc_actual.CVA_vars_girf.irf.(f); 
    S1_60_pc_actual.aLP_MSPE_CVA_girf.mu.(f) = mean(S1_60_pc_actual.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_actual.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S1_60_pc_actual.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S1_60_pc_actual.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S1_60_pc_actual.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S1_60_pc_actual','-append')
%% ------------------------------------------------------- 
% 1-3. the rolling window (30 months) correlation between inflation (yoy) and
% consumption growth (actual)
% S1_30_id_actual
% -------------------------------------------------------
    
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\actual\lag12\Bloom_S1_30_id_actual.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S1_30_id_actual.(names{i}).irf.(f) = irfpaths.(names{i}).S1_30_id.(f);
        S1_30_id_actual.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S1_30_id.(f),2);
        S1_30_id_actual.(names{i}).band_hi.(f) = band_hi.(names{i}).S1_30_id.(f);
        S1_30_id_actual.(names{i}).band_lo.(f) = band_lo.(names{i}).S1_30_id.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S1_30_id.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S1_30_id_actual.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_30_id.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S1_30_id_actual.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_30_id.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S1_30_id_actual.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S1_30_id.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S1_30_id_actual.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S1_30_id.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S1_30_id_actual.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S1_30_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S1_30_id_actual.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_actual.GMA_lps.irf.(f) = wirf_val;% H x B
    S1_30_id_actual.GMA_lps.mu.(f) = mean(S1_30_id_actual.GMA_lps.irf.(f), 2, "omitmissing");
    S1_30_id_actual.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_actual.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S1_30_id_actual.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_actual.CVA_lps.irf.(f) = wirf_val;% H x B
    S1_30_id_actual.CVA_lps.mu.(f) = mean(S1_30_id_actual.CVA_lps.irf.(f), 2, "omitmissing");
    S1_30_id_actual.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_actual.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S1_30_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_30_id_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_actual.GMA_vars.irf.(f) = wirf_val;% H x B
    S1_30_id_actual.GMA_vars.mu.(f) = mean(S1_30_id_actual.GMA_vars.irf.(f), 2, "omitmissing");
    S1_30_id_actual.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_actual.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_30_id_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_actual.CVA_vars.irf.(f) = wirf_val;
    S1_30_id_actual.CVA_vars.mu.(f) = mean(S1_30_id_actual.CVA_vars.irf.(f), 2, "omitmissing");
    S1_30_id_actual.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_actual.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S1_30_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_30_id_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_actual.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S1_30_id_actual.GMA_vars_girf.mu.(f) = mean(S1_30_id_actual.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S1_30_id_actual.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_actual.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_30_id_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_actual.CVA_vars_girf.irf.(f) = wirf_val;
    S1_30_id_actual.CVA_vars_girf.mu.(f) = mean(S1_30_id_actual.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S1_30_id_actual.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_actual.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S1_30_id.state0_vol_empl = tss.S1_30_id.state0_lps_empl;
    tss.S1_30_id.state0_ff_empl = tss.S1_30_id.state0_lps_empl;
    tss.S1_30_id.state1_vol_empl = tss.S1_30_id.state1_lps_empl;
    tss.S1_30_id.state1_ff_empl = tss.S1_30_id.state1_lps_empl;
    tss.S1_30_id.state0_vol_ip = tss.S1_30_id.state0_lps_ip;
    tss.S1_30_id.state0_ff_ip = tss.S1_30_id.state0_lps_ip;
    tss.S1_30_id.state1_vol_ip = tss.S1_30_id.state1_lps_ip;
    tss.S1_30_id.state1_ff_ip = tss.S1_30_id.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S1_30_id.(f) ./ tss.S1_30_id.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S1_30_id.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S1_30_id.state0_vol_empl = tss.S1_30_id.state0_vars_empl;
    tss.S1_30_id.state0_ff_empl = tss.S1_30_id.state0_vars_empl;
    tss.S1_30_id.state1_vol_empl = tss.S1_30_id.state1_vars_empl;
    tss.S1_30_id.state1_ff_empl = tss.S1_30_id.state1_vars_empl;
    tss.S1_30_id.state0_vol_ip = tss.S1_30_id.state0_vars_ip;
    tss.S1_30_id.state0_ff_ip = tss.S1_30_id.state0_vars_ip;
    tss.S1_30_id.state1_vol_ip = tss.S1_30_id.state1_vars_ip;
    tss.S1_30_id.state1_ff_ip = tss.S1_30_id.state1_vars_ip;
    tss.S1_30_id.linear_vol_empl = tss.S1_30_id.linear_vars_empl;
    tss.S1_30_id.linear_vol_ip = tss.S1_30_id.linear_vars_ip;
    tss.S1_30_id.linear_ff_empl = tss.S1_30_id.linear_vars_empl;
    tss.S1_30_id.linear_ff_ip = tss.S1_30_id.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S1_30_id.(f) ./ tss.S1_30_id.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S1_30_id.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S1_30_id_actual.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S1_30_id_actual.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S1_30_id_actual.aLP_RS.alpha.(f) .* S1_30_id_actual.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S1_30_id_actual.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_30_id_actual.GMA_vars.weight.(f)'; % H x 2
    S1_30_id_actual.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S1_30_id_actual.aLP_RS.alpha.(f) .* S1_30_id_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_30_id_actual.CVA_vars.weight.(f)';
    S1_30_id_actual.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S1_30_id_actual.aLP_MSPE.alpha.(f) .* S1_30_id_actual.GMA_lps.weight.(f)';
    alpha_var = 1 - S1_30_id_actual.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_30_id_actual.GMA_vars.weight.(f)';
    S1_30_id_actual.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S1_30_id_actual.aLP_MSPE.alpha.(f) .* S1_30_id_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_30_id_actual.CVA_vars.weight.(f)';
    S1_30_id_actual.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S1_30_id_actual.aLP_RS_GMA.irf.(f) = S1_30_id_actual.aLP_RS.alpha.(f) .* S1_30_id_actual.GMA_lps.irf.(f)...
    + (1-S1_30_id_actual.aLP_RS.alpha.(f)) .* S1_30_id_actual.GMA_vars.irf.(f); 
    S1_30_id_actual.aLP_RS_GMA.mu.(f) = mean(S1_30_id_actual.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S1_30_id_actual.aLP_RS_GMA.band_lo.(f) = prctile(S1_30_id_actual.aLP_RS_GMA.irf.(f),pct_lo,2);
    S1_30_id_actual.aLP_RS_GMA.band_hi.(f) = prctile(S1_30_id_actual.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S1_30_id_actual.aLP_RS_CVA.irf.(f) = S1_30_id_actual.aLP_RS.alpha.(f) .* S1_30_id_actual.CVA_lps.irf.(f)...
    + (1-S1_30_id_actual.aLP_RS.alpha.(f)) .* S1_30_id_actual.CVA_vars.irf.(f); 
    S1_30_id_actual.aLP_RS_CVA.mu.(f) = mean(S1_30_id_actual.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S1_30_id_actual.aLP_RS_CVA.band_lo.(f) = prctile(S1_30_id_actual.aLP_RS_CVA.irf.(f),pct_lo,2);
    S1_30_id_actual.aLP_RS_CVA.band_hi.(f) = prctile(S1_30_id_actual.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S1_30_id_actual.aLP_MSPE_GMA.irf.(f) = S1_30_id_actual.aLP_MSPE.alpha.(f) .* S1_30_id_actual.GMA_lps.irf.(f)...
    + (1-S1_30_id_actual.aLP_MSPE.alpha.(f)) .* S1_30_id_actual.GMA_vars.irf.(f); 
    S1_30_id_actual.aLP_MSPE_GMA.mu.(f) = mean(S1_30_id_actual.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S1_30_id_actual.aLP_MSPE_GMA.band_lo.(f) = prctile(S1_30_id_actual.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S1_30_id_actual.aLP_MSPE_GMA.band_hi.(f) = prctile(S1_30_id_actual.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S1_30_id_actual.aLP_MSPE_CVA.irf.(f) = S1_30_id_actual.aLP_MSPE.alpha.(f) .* S1_30_id_actual.CVA_lps.irf.(f)...
    + (1-S1_30_id_actual.aLP_MSPE.alpha.(f)) .* S1_30_id_actual.CVA_vars.irf.(f); 
    S1_30_id_actual.aLP_MSPE_CVA.mu.(f) = mean(S1_30_id_actual.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S1_30_id_actual.aLP_MSPE_CVA.band_lo.(f) = prctile(S1_30_id_actual.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S1_30_id_actual.aLP_MSPE_CVA.band_hi.(f) = prctile(S1_30_id_actual.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S1_30_id_actual.aLP_RS_GMA_girf.irf.(f) = S1_30_id_actual.aLP_RS.alpha.(f) .* S1_30_id_actual.GMA_lps.irf.(f)...
    + (1-S1_30_id_actual.aLP_RS.alpha.(f)) .* S1_30_id_actual.GMA_vars_girf.irf.(f); 
    S1_30_id_actual.aLP_RS_GMA_girf.mu.(f) = mean(S1_30_id_actual.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S1_30_id_actual.aLP_RS_GMA_girf.band_lo.(f) = prctile(S1_30_id_actual.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S1_30_id_actual.aLP_RS_GMA_girf.band_hi.(f) = prctile(S1_30_id_actual.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_30_id_actual.aLP_RS_CVA_girf.irf.(f) = S1_30_id_actual.aLP_RS.alpha.(f) .* S1_30_id_actual.CVA_lps.irf.(f)...
    + (1-S1_30_id_actual.aLP_RS.alpha.(f)) .* S1_30_id_actual.CVA_vars_girf.irf.(f); 
    S1_30_id_actual.aLP_RS_CVA_girf.mu.(f) = mean(S1_30_id_actual.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S1_30_id_actual.aLP_RS_CVA_girf.band_lo.(f) = prctile(S1_30_id_actual.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S1_30_id_actual.aLP_RS_CVA_girf.band_hi.(f) = prctile(S1_30_id_actual.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S1_30_id_actual.aLP_MSPE_GMA_girf.irf.(f) = S1_30_id_actual.aLP_MSPE.alpha.(f) .* S1_30_id_actual.GMA_lps.irf.(f)...
    + (1-S1_30_id_actual.aLP_MSPE.alpha.(f)) .* S1_30_id_actual.GMA_vars_girf.irf.(f); 
    S1_30_id_actual.aLP_MSPE_GMA_girf.mu.(f) = mean(S1_30_id_actual.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S1_30_id_actual.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S1_30_id_actual.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S1_30_id_actual.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S1_30_id_actual.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_30_id_actual.aLP_MSPE_CVA_girf.irf.(f) = S1_30_id_actual.aLP_MSPE.alpha.(f) .* S1_30_id_actual.CVA_lps.irf.(f)...
    + (1-S1_30_id_actual.aLP_MSPE.alpha.(f)) .* S1_30_id_actual.CVA_vars_girf.irf.(f); 
    S1_30_id_actual.aLP_MSPE_CVA_girf.mu.(f) = mean(S1_30_id_actual.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S1_30_id_actual.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S1_30_id_actual.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S1_30_id_actual.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S1_30_id_actual.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S1_30_id_actual','-append')
%% ------------------------------------------------------- 
% 1-4. the rolling window (30 months) correlation between inflation (yoy) and
% consumption growth (per capita)
% S1_30_pc_actual
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\actual\lag12\Bloom_S1_30_pc_actual.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S1_30_pc_actual.(names{i}).irf.(f) = irfpaths.(names{i}).S1_30_pc.(f);
        S1_30_pc_actual.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S1_30_pc.(f),2);
        S1_30_pc_actual.(names{i}).band_hi.(f) = band_hi.(names{i}).S1_30_pc.(f);
        S1_30_pc_actual.(names{i}).band_lo.(f) = band_lo.(names{i}).S1_30_pc.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S1_30_pc.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S1_30_pc_actual.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_30_pc.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S1_30_pc_actual.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_30_pc.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S1_30_pc_actual.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S1_30_pc.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S1_30_pc_actual.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S1_30_pc.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S1_30_pc_actual.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S1_30_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S1_30_pc_actual.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_actual.GMA_lps.irf.(f) = wirf_val;% H x B
    S1_30_pc_actual.GMA_lps.mu.(f) = mean(S1_30_pc_actual.GMA_lps.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_actual.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S1_30_pc_actual.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_actual.CVA_lps.irf.(f) = wirf_val;% H x B
    S1_30_pc_actual.CVA_lps.mu.(f) = mean(S1_30_pc_actual.CVA_lps.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_actual.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S1_30_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_30_pc_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_actual.GMA_vars.irf.(f) = wirf_val;% H x B
    S1_30_pc_actual.GMA_vars.mu.(f) = mean(S1_30_pc_actual.GMA_vars.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_actual.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_30_pc_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_actual.CVA_vars.irf.(f) = wirf_val;
    S1_30_pc_actual.CVA_vars.mu.(f) = mean(S1_30_pc_actual.CVA_vars.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_actual.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S1_30_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_30_pc_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_actual.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S1_30_pc_actual.GMA_vars_girf.mu.(f) = mean(S1_30_pc_actual.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_actual.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_30_pc_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_actual.CVA_vars_girf.irf.(f) = wirf_val;
    S1_30_pc_actual.CVA_vars_girf.mu.(f) = mean(S1_30_pc_actual.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_actual.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S1_30_pc.state0_vol_empl = tss.S1_30_pc.state0_lps_empl;
    tss.S1_30_pc.state0_ff_empl = tss.S1_30_pc.state0_lps_empl;
    tss.S1_30_pc.state1_vol_empl = tss.S1_30_pc.state1_lps_empl;
    tss.S1_30_pc.state1_ff_empl = tss.S1_30_pc.state1_lps_empl;
    tss.S1_30_pc.state0_vol_ip = tss.S1_30_pc.state0_lps_ip;
    tss.S1_30_pc.state0_ff_ip = tss.S1_30_pc.state0_lps_ip;
    tss.S1_30_pc.state1_vol_ip = tss.S1_30_pc.state1_lps_ip;
    tss.S1_30_pc.state1_ff_ip = tss.S1_30_pc.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S1_30_pc.(f) ./ tss.S1_30_pc.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S1_30_pc.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S1_30_pc.state0_vol_empl = tss.S1_30_pc.state0_vars_empl;
    tss.S1_30_pc.state0_ff_empl = tss.S1_30_pc.state0_vars_empl;
    tss.S1_30_pc.state1_vol_empl = tss.S1_30_pc.state1_vars_empl;
    tss.S1_30_pc.state1_ff_empl = tss.S1_30_pc.state1_vars_empl;
    tss.S1_30_pc.state0_vol_ip = tss.S1_30_pc.state0_vars_ip;
    tss.S1_30_pc.state0_ff_ip = tss.S1_30_pc.state0_vars_ip;
    tss.S1_30_pc.state1_vol_ip = tss.S1_30_pc.state1_vars_ip;
    tss.S1_30_pc.state1_ff_ip = tss.S1_30_pc.state1_vars_ip;
    tss.S1_30_pc.linear_vol_empl = tss.S1_30_pc.linear_vars_empl;
    tss.S1_30_pc.linear_vol_ip = tss.S1_30_pc.linear_vars_ip;
    tss.S1_30_pc.linear_ff_empl = tss.S1_30_pc.linear_vars_empl;
    tss.S1_30_pc.linear_ff_ip = tss.S1_30_pc.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S1_30_pc.(f) ./ tss.S1_30_pc.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S1_30_pc.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S1_30_pc_actual.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S1_30_pc_actual.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S1_30_pc_actual.aLP_RS.alpha.(f) .* S1_30_pc_actual.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S1_30_pc_actual.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_30_pc_actual.GMA_vars.weight.(f)'; % H x 2
    S1_30_pc_actual.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S1_30_pc_actual.aLP_RS.alpha.(f) .* S1_30_pc_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_30_pc_actual.CVA_vars.weight.(f)';
    S1_30_pc_actual.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S1_30_pc_actual.aLP_MSPE.alpha.(f) .* S1_30_pc_actual.GMA_lps.weight.(f)';
    alpha_var = 1 - S1_30_pc_actual.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_30_pc_actual.GMA_vars.weight.(f)';
    S1_30_pc_actual.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S1_30_pc_actual.aLP_MSPE.alpha.(f) .* S1_30_pc_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_30_pc_actual.CVA_vars.weight.(f)';
    S1_30_pc_actual.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S1_30_pc_actual.aLP_RS_GMA.irf.(f) = S1_30_pc_actual.aLP_RS.alpha.(f) .* S1_30_pc_actual.GMA_lps.irf.(f)...
    + (1-S1_30_pc_actual.aLP_RS.alpha.(f)) .* S1_30_pc_actual.GMA_vars.irf.(f); 
    S1_30_pc_actual.aLP_RS_GMA.mu.(f) = mean(S1_30_pc_actual.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.aLP_RS_GMA.band_lo.(f) = prctile(S1_30_pc_actual.aLP_RS_GMA.irf.(f),pct_lo,2);
    S1_30_pc_actual.aLP_RS_GMA.band_hi.(f) = prctile(S1_30_pc_actual.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S1_30_pc_actual.aLP_RS_CVA.irf.(f) = S1_30_pc_actual.aLP_RS.alpha.(f) .* S1_30_pc_actual.CVA_lps.irf.(f)...
    + (1-S1_30_pc_actual.aLP_RS.alpha.(f)) .* S1_30_pc_actual.CVA_vars.irf.(f); 
    S1_30_pc_actual.aLP_RS_CVA.mu.(f) = mean(S1_30_pc_actual.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.aLP_RS_CVA.band_lo.(f) = prctile(S1_30_pc_actual.aLP_RS_CVA.irf.(f),pct_lo,2);
    S1_30_pc_actual.aLP_RS_CVA.band_hi.(f) = prctile(S1_30_pc_actual.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S1_30_pc_actual.aLP_MSPE_GMA.irf.(f) = S1_30_pc_actual.aLP_MSPE.alpha.(f) .* S1_30_pc_actual.GMA_lps.irf.(f)...
    + (1-S1_30_pc_actual.aLP_MSPE.alpha.(f)) .* S1_30_pc_actual.GMA_vars.irf.(f); 
    S1_30_pc_actual.aLP_MSPE_GMA.mu.(f) = mean(S1_30_pc_actual.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.aLP_MSPE_GMA.band_lo.(f) = prctile(S1_30_pc_actual.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S1_30_pc_actual.aLP_MSPE_GMA.band_hi.(f) = prctile(S1_30_pc_actual.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S1_30_pc_actual.aLP_MSPE_CVA.irf.(f) = S1_30_pc_actual.aLP_MSPE.alpha.(f) .* S1_30_pc_actual.CVA_lps.irf.(f)...
    + (1-S1_30_pc_actual.aLP_MSPE.alpha.(f)) .* S1_30_pc_actual.CVA_vars.irf.(f); 
    S1_30_pc_actual.aLP_MSPE_CVA.mu.(f) = mean(S1_30_pc_actual.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.aLP_MSPE_CVA.band_lo.(f) = prctile(S1_30_pc_actual.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S1_30_pc_actual.aLP_MSPE_CVA.band_hi.(f) = prctile(S1_30_pc_actual.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S1_30_pc_actual.aLP_RS_GMA_girf.irf.(f) = S1_30_pc_actual.aLP_RS.alpha.(f) .* S1_30_pc_actual.GMA_lps.irf.(f)...
    + (1-S1_30_pc_actual.aLP_RS.alpha.(f)) .* S1_30_pc_actual.GMA_vars_girf.irf.(f); 
    S1_30_pc_actual.aLP_RS_GMA_girf.mu.(f) = mean(S1_30_pc_actual.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.aLP_RS_GMA_girf.band_lo.(f) = prctile(S1_30_pc_actual.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S1_30_pc_actual.aLP_RS_GMA_girf.band_hi.(f) = prctile(S1_30_pc_actual.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_30_pc_actual.aLP_RS_CVA_girf.irf.(f) = S1_30_pc_actual.aLP_RS.alpha.(f) .* S1_30_pc_actual.CVA_lps.irf.(f)...
    + (1-S1_30_pc_actual.aLP_RS.alpha.(f)) .* S1_30_pc_actual.CVA_vars_girf.irf.(f); 
    S1_30_pc_actual.aLP_RS_CVA_girf.mu.(f) = mean(S1_30_pc_actual.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.aLP_RS_CVA_girf.band_lo.(f) = prctile(S1_30_pc_actual.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S1_30_pc_actual.aLP_RS_CVA_girf.band_hi.(f) = prctile(S1_30_pc_actual.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S1_30_pc_actual.aLP_MSPE_GMA_girf.irf.(f) = S1_30_pc_actual.aLP_MSPE.alpha.(f) .* S1_30_pc_actual.GMA_lps.irf.(f)...
    + (1-S1_30_pc_actual.aLP_MSPE.alpha.(f)) .* S1_30_pc_actual.GMA_vars_girf.irf.(f); 
    S1_30_pc_actual.aLP_MSPE_GMA_girf.mu.(f) = mean(S1_30_pc_actual.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S1_30_pc_actual.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S1_30_pc_actual.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S1_30_pc_actual.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_30_pc_actual.aLP_MSPE_CVA_girf.irf.(f) = S1_30_pc_actual.aLP_MSPE.alpha.(f) .* S1_30_pc_actual.CVA_lps.irf.(f)...
    + (1-S1_30_pc_actual.aLP_MSPE.alpha.(f)) .* S1_30_pc_actual.CVA_vars_girf.irf.(f); 
    S1_30_pc_actual.aLP_MSPE_CVA_girf.mu.(f) = mean(S1_30_pc_actual.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_actual.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S1_30_pc_actual.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S1_30_pc_actual.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S1_30_pc_actual.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S1_30_pc_actual','-append')
%% -------------------------------------------------------
% 2. Anchored vs. Unanchored expectations
% 2-1. Michigan Household Survey IQR from jan 1987
% S2_78_hh_actual
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\actual\lag12\Bloom_S2_78_hh_actual.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S2_78_hh_actual.(names{i}).irf.(f) = irfpaths.(names{i}).S2_78_hh.(f);
        S2_78_hh_actual.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S2_78_hh.(f),2);
        S2_78_hh_actual.(names{i}).band_hi.(f) = band_hi.(names{i}).S2_78_hh.(f);
        S2_78_hh_actual.(names{i}).band_lo.(f) = band_lo.(names{i}).S2_78_hh.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S2_78_hh.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S2_78_hh_actual.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_78_hh.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S2_78_hh_actual.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_78_hh.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S2_78_hh_actual.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S2_78_hh.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S2_78_hh_actual.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S2_78_hh.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S2_78_hh_actual.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S2_78_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S2_78_hh_actual.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_actual.GMA_lps.irf.(f) = wirf_val;% H x B
    S2_78_hh_actual.GMA_lps.mu.(f) = mean(S2_78_hh_actual.GMA_lps.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_actual.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S2_78_hh_actual.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_actual.CVA_lps.irf.(f) = wirf_val;% H x B
    S2_78_hh_actual.CVA_lps.mu.(f) = mean(S2_78_hh_actual.CVA_lps.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_actual.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S2_78_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_78_hh_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_actual.GMA_vars.irf.(f) = wirf_val;% H x B
    S2_78_hh_actual.GMA_vars.mu.(f) = mean(S2_78_hh_actual.GMA_vars.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_actual.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_78_hh_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_actual.CVA_vars.irf.(f) = wirf_val;
    S2_78_hh_actual.CVA_vars.mu.(f) = mean(S2_78_hh_actual.CVA_vars.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_actual.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S2_78_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_78_hh_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_actual.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S2_78_hh_actual.GMA_vars_girf.mu.(f) = mean(S2_78_hh_actual.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_actual.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_78_hh_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_actual.CVA_vars_girf.irf.(f) = wirf_val;
    S2_78_hh_actual.CVA_vars_girf.mu.(f) = mean(S2_78_hh_actual.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_actual.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S2_78_hh.state0_vol_empl = tss.S2_78_hh.state0_lps_empl;
    tss.S2_78_hh.state0_ff_empl = tss.S2_78_hh.state0_lps_empl;
    tss.S2_78_hh.state1_vol_empl = tss.S2_78_hh.state1_lps_empl;
    tss.S2_78_hh.state1_ff_empl = tss.S2_78_hh.state1_lps_empl;
    tss.S2_78_hh.state0_vol_ip = tss.S2_78_hh.state0_lps_ip;
    tss.S2_78_hh.state0_ff_ip = tss.S2_78_hh.state0_lps_ip;
    tss.S2_78_hh.state1_vol_ip = tss.S2_78_hh.state1_lps_ip;
    tss.S2_78_hh.state1_ff_ip = tss.S2_78_hh.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S2_78_hh.(f) ./ tss.S2_78_hh.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S2_78_hh.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S2_78_hh.state0_vol_empl = tss.S2_78_hh.state0_vars_empl;
    tss.S2_78_hh.state0_ff_empl = tss.S2_78_hh.state0_vars_empl;
    tss.S2_78_hh.state1_vol_empl = tss.S2_78_hh.state1_vars_empl;
    tss.S2_78_hh.state1_ff_empl = tss.S2_78_hh.state1_vars_empl;
    tss.S2_78_hh.state0_vol_ip = tss.S2_78_hh.state0_vars_ip;
    tss.S2_78_hh.state0_ff_ip = tss.S2_78_hh.state0_vars_ip;
    tss.S2_78_hh.state1_vol_ip = tss.S2_78_hh.state1_vars_ip;
    tss.S2_78_hh.state1_ff_ip = tss.S2_78_hh.state1_vars_ip;
    tss.S2_78_hh.linear_vol_empl = tss.S2_78_hh.linear_vars_empl;
    tss.S2_78_hh.linear_vol_ip = tss.S2_78_hh.linear_vars_ip;
    tss.S2_78_hh.linear_ff_empl = tss.S2_78_hh.linear_vars_empl;
    tss.S2_78_hh.linear_ff_ip = tss.S2_78_hh.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S2_78_hh.(f) ./ tss.S2_78_hh.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S2_78_hh.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S2_78_hh_actual.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S2_78_hh_actual.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S2_78_hh_actual.aLP_RS.alpha.(f) .* S2_78_hh_actual.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S2_78_hh_actual.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_78_hh_actual.GMA_vars.weight.(f)'; % H x 2
    S2_78_hh_actual.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S2_78_hh_actual.aLP_RS.alpha.(f) .* S2_78_hh_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_78_hh_actual.CVA_vars.weight.(f)';
    S2_78_hh_actual.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S2_78_hh_actual.aLP_MSPE.alpha.(f) .* S2_78_hh_actual.GMA_lps.weight.(f)';
    alpha_var = 1 - S2_78_hh_actual.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_78_hh_actual.GMA_vars.weight.(f)';
    S2_78_hh_actual.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S2_78_hh_actual.aLP_MSPE.alpha.(f) .* S2_78_hh_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_78_hh_actual.CVA_vars.weight.(f)';
    S2_78_hh_actual.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S2_78_hh_actual.aLP_RS_GMA.irf.(f) = S2_78_hh_actual.aLP_RS.alpha.(f) .* S2_78_hh_actual.GMA_lps.irf.(f)...
    + (1-S2_78_hh_actual.aLP_RS.alpha.(f)) .* S2_78_hh_actual.GMA_vars.irf.(f); 
    S2_78_hh_actual.aLP_RS_GMA.mu.(f) = mean(S2_78_hh_actual.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.aLP_RS_GMA.band_lo.(f) = prctile(S2_78_hh_actual.aLP_RS_GMA.irf.(f),pct_lo,2);
    S2_78_hh_actual.aLP_RS_GMA.band_hi.(f) = prctile(S2_78_hh_actual.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S2_78_hh_actual.aLP_RS_CVA.irf.(f) = S2_78_hh_actual.aLP_RS.alpha.(f) .* S2_78_hh_actual.CVA_lps.irf.(f)...
    + (1-S2_78_hh_actual.aLP_RS.alpha.(f)) .* S2_78_hh_actual.CVA_vars.irf.(f); 
    S2_78_hh_actual.aLP_RS_CVA.mu.(f) = mean(S2_78_hh_actual.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.aLP_RS_CVA.band_lo.(f) = prctile(S2_78_hh_actual.aLP_RS_CVA.irf.(f),pct_lo,2);
    S2_78_hh_actual.aLP_RS_CVA.band_hi.(f) = prctile(S2_78_hh_actual.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S2_78_hh_actual.aLP_MSPE_GMA.irf.(f) = S2_78_hh_actual.aLP_MSPE.alpha.(f) .* S2_78_hh_actual.GMA_lps.irf.(f)...
    + (1-S2_78_hh_actual.aLP_MSPE.alpha.(f)) .* S2_78_hh_actual.GMA_vars.irf.(f); 
    S2_78_hh_actual.aLP_MSPE_GMA.mu.(f) = mean(S2_78_hh_actual.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.aLP_MSPE_GMA.band_lo.(f) = prctile(S2_78_hh_actual.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S2_78_hh_actual.aLP_MSPE_GMA.band_hi.(f) = prctile(S2_78_hh_actual.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S2_78_hh_actual.aLP_MSPE_CVA.irf.(f) = S2_78_hh_actual.aLP_MSPE.alpha.(f) .* S2_78_hh_actual.CVA_lps.irf.(f)...
    + (1-S2_78_hh_actual.aLP_MSPE.alpha.(f)) .* S2_78_hh_actual.CVA_vars.irf.(f); 
    S2_78_hh_actual.aLP_MSPE_CVA.mu.(f) = mean(S2_78_hh_actual.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.aLP_MSPE_CVA.band_lo.(f) = prctile(S2_78_hh_actual.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S2_78_hh_actual.aLP_MSPE_CVA.band_hi.(f) = prctile(S2_78_hh_actual.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S2_78_hh_actual.aLP_RS_GMA_girf.irf.(f) = S2_78_hh_actual.aLP_RS.alpha.(f) .* S2_78_hh_actual.GMA_lps.irf.(f)...
    + (1-S2_78_hh_actual.aLP_RS.alpha.(f)) .* S2_78_hh_actual.GMA_vars_girf.irf.(f); 
    S2_78_hh_actual.aLP_RS_GMA_girf.mu.(f) = mean(S2_78_hh_actual.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.aLP_RS_GMA_girf.band_lo.(f) = prctile(S2_78_hh_actual.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S2_78_hh_actual.aLP_RS_GMA_girf.band_hi.(f) = prctile(S2_78_hh_actual.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_78_hh_actual.aLP_RS_CVA_girf.irf.(f) = S2_78_hh_actual.aLP_RS.alpha.(f) .* S2_78_hh_actual.CVA_lps.irf.(f)...
    + (1-S2_78_hh_actual.aLP_RS.alpha.(f)) .* S2_78_hh_actual.CVA_vars_girf.irf.(f); 
    S2_78_hh_actual.aLP_RS_CVA_girf.mu.(f) = mean(S2_78_hh_actual.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.aLP_RS_CVA_girf.band_lo.(f) = prctile(S2_78_hh_actual.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S2_78_hh_actual.aLP_RS_CVA_girf.band_hi.(f) = prctile(S2_78_hh_actual.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S2_78_hh_actual.aLP_MSPE_GMA_girf.irf.(f) = S2_78_hh_actual.aLP_MSPE.alpha.(f) .* S2_78_hh_actual.GMA_lps.irf.(f)...
    + (1-S2_78_hh_actual.aLP_MSPE.alpha.(f)) .* S2_78_hh_actual.GMA_vars_girf.irf.(f); 
    S2_78_hh_actual.aLP_MSPE_GMA_girf.mu.(f) = mean(S2_78_hh_actual.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S2_78_hh_actual.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S2_78_hh_actual.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S2_78_hh_actual.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_78_hh_actual.aLP_MSPE_CVA_girf.irf.(f) = S2_78_hh_actual.aLP_MSPE.alpha.(f) .* S2_78_hh_actual.CVA_lps.irf.(f)...
    + (1-S2_78_hh_actual.aLP_MSPE.alpha.(f)) .* S2_78_hh_actual.CVA_vars_girf.irf.(f); 
    S2_78_hh_actual.aLP_MSPE_CVA_girf.mu.(f) = mean(S2_78_hh_actual.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_actual.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S2_78_hh_actual.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S2_78_hh_actual.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S2_78_hh_actual.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end
%%
save("Combine_Results_Emp.mat",'S2_78_hh_actual','-append')
%% ------------------------------------------------------
% 2-2. Michigan Household Survey IQR from 1981
% S2_81_hh_actual
% -------------------------------------------------------   
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\actual\lag12\Bloom_S2_81_hh_actual.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S2_81_hh_actual.(names{i}).irf.(f) = irfpaths.(names{i}).S2_81_hh.(f);
        S2_81_hh_actual.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S2_81_hh.(f),2);
        S2_81_hh_actual.(names{i}).band_hi.(f) = band_hi.(names{i}).S2_81_hh.(f);
        S2_81_hh_actual.(names{i}).band_lo.(f) = band_lo.(names{i}).S2_81_hh.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S2_81_hh.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S2_81_hh_actual.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_81_hh.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S2_81_hh_actual.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_81_hh.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S2_81_hh_actual.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S2_81_hh.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S2_81_hh_actual.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S2_81_hh.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S2_81_hh_actual.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S2_81_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S2_81_hh_actual.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_actual.GMA_lps.irf.(f) = wirf_val;% H x B
    S2_81_hh_actual.GMA_lps.mu.(f) = mean(S2_81_hh_actual.GMA_lps.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_actual.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S2_81_hh_actual.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_actual.CVA_lps.irf.(f) = wirf_val;% H x B
    S2_81_hh_actual.CVA_lps.mu.(f) = mean(S2_81_hh_actual.CVA_lps.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_actual.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S2_81_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_81_hh_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_actual.GMA_vars.irf.(f) = wirf_val;% H x B
    S2_81_hh_actual.GMA_vars.mu.(f) = mean(S2_81_hh_actual.GMA_vars.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_actual.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_81_hh_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_actual.CVA_vars.irf.(f) = wirf_val;
    S2_81_hh_actual.CVA_vars.mu.(f) = mean(S2_81_hh_actual.CVA_vars.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_actual.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S2_81_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_81_hh_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_actual.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S2_81_hh_actual.GMA_vars_girf.mu.(f) = mean(S2_81_hh_actual.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_actual.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_81_hh_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_actual.CVA_vars_girf.irf.(f) = wirf_val;
    S2_81_hh_actual.CVA_vars_girf.mu.(f) = mean(S2_81_hh_actual.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_actual.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S2_81_hh.state0_vol_empl = tss.S2_81_hh.state0_lps_empl;
    tss.S2_81_hh.state0_ff_empl = tss.S2_81_hh.state0_lps_empl;
    tss.S2_81_hh.state1_vol_empl = tss.S2_81_hh.state1_lps_empl;
    tss.S2_81_hh.state1_ff_empl = tss.S2_81_hh.state1_lps_empl;
    tss.S2_81_hh.state0_vol_ip = tss.S2_81_hh.state0_lps_ip;
    tss.S2_81_hh.state0_ff_ip = tss.S2_81_hh.state0_lps_ip;
    tss.S2_81_hh.state1_vol_ip = tss.S2_81_hh.state1_lps_ip;
    tss.S2_81_hh.state1_ff_ip = tss.S2_81_hh.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S2_81_hh.(f) ./ tss.S2_81_hh.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S2_81_hh.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S2_81_hh.state0_vol_empl = tss.S2_81_hh.state0_vars_empl;
    tss.S2_81_hh.state0_ff_empl = tss.S2_81_hh.state0_vars_empl;
    tss.S2_81_hh.state1_vol_empl = tss.S2_81_hh.state1_vars_empl;
    tss.S2_81_hh.state1_ff_empl = tss.S2_81_hh.state1_vars_empl;
    tss.S2_81_hh.state0_vol_ip = tss.S2_81_hh.state0_vars_ip;
    tss.S2_81_hh.state0_ff_ip = tss.S2_81_hh.state0_vars_ip;
    tss.S2_81_hh.state1_vol_ip = tss.S2_81_hh.state1_vars_ip;
    tss.S2_81_hh.state1_ff_ip = tss.S2_81_hh.state1_vars_ip;
    tss.S2_81_hh.linear_vol_empl = tss.S2_81_hh.linear_vars_empl;
    tss.S2_81_hh.linear_vol_ip = tss.S2_81_hh.linear_vars_ip;
    tss.S2_81_hh.linear_ff_empl = tss.S2_81_hh.linear_vars_empl;
    tss.S2_81_hh.linear_ff_ip = tss.S2_81_hh.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S2_81_hh.(f) ./ tss.S2_81_hh.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S2_81_hh.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S2_81_hh_actual.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S2_81_hh_actual.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S2_81_hh_actual.aLP_RS.alpha.(f) .* S2_81_hh_actual.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S2_81_hh_actual.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_81_hh_actual.GMA_vars.weight.(f)'; % H x 2
    S2_81_hh_actual.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S2_81_hh_actual.aLP_RS.alpha.(f) .* S2_81_hh_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_81_hh_actual.CVA_vars.weight.(f)';
    S2_81_hh_actual.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S2_81_hh_actual.aLP_MSPE.alpha.(f) .* S2_81_hh_actual.GMA_lps.weight.(f)';
    alpha_var = 1 - S2_81_hh_actual.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_81_hh_actual.GMA_vars.weight.(f)';
    S2_81_hh_actual.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S2_81_hh_actual.aLP_MSPE.alpha.(f) .* S2_81_hh_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_81_hh_actual.CVA_vars.weight.(f)';
    S2_81_hh_actual.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S2_81_hh_actual.aLP_RS_GMA.irf.(f) = S2_81_hh_actual.aLP_RS.alpha.(f) .* S2_81_hh_actual.GMA_lps.irf.(f)...
    + (1-S2_81_hh_actual.aLP_RS.alpha.(f)) .* S2_81_hh_actual.GMA_vars.irf.(f); 
    S2_81_hh_actual.aLP_RS_GMA.mu.(f) = mean(S2_81_hh_actual.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.aLP_RS_GMA.band_lo.(f) = prctile(S2_81_hh_actual.aLP_RS_GMA.irf.(f),pct_lo,2);
    S2_81_hh_actual.aLP_RS_GMA.band_hi.(f) = prctile(S2_81_hh_actual.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S2_81_hh_actual.aLP_RS_CVA.irf.(f) = S2_81_hh_actual.aLP_RS.alpha.(f) .* S2_81_hh_actual.CVA_lps.irf.(f)...
    + (1-S2_81_hh_actual.aLP_RS.alpha.(f)) .* S2_81_hh_actual.CVA_vars.irf.(f); 
    S2_81_hh_actual.aLP_RS_CVA.mu.(f) = mean(S2_81_hh_actual.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.aLP_RS_CVA.band_lo.(f) = prctile(S2_81_hh_actual.aLP_RS_CVA.irf.(f),pct_lo,2);
    S2_81_hh_actual.aLP_RS_CVA.band_hi.(f) = prctile(S2_81_hh_actual.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S2_81_hh_actual.aLP_MSPE_GMA.irf.(f) = S2_81_hh_actual.aLP_MSPE.alpha.(f) .* S2_81_hh_actual.GMA_lps.irf.(f)...
    + (1-S2_81_hh_actual.aLP_MSPE.alpha.(f)) .* S2_81_hh_actual.GMA_vars.irf.(f); 
    S2_81_hh_actual.aLP_MSPE_GMA.mu.(f) = mean(S2_81_hh_actual.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.aLP_MSPE_GMA.band_lo.(f) = prctile(S2_81_hh_actual.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S2_81_hh_actual.aLP_MSPE_GMA.band_hi.(f) = prctile(S2_81_hh_actual.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S2_81_hh_actual.aLP_MSPE_CVA.irf.(f) = S2_81_hh_actual.aLP_MSPE.alpha.(f) .* S2_81_hh_actual.CVA_lps.irf.(f)...
    + (1-S2_81_hh_actual.aLP_MSPE.alpha.(f)) .* S2_81_hh_actual.CVA_vars.irf.(f); 
    S2_81_hh_actual.aLP_MSPE_CVA.mu.(f) = mean(S2_81_hh_actual.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.aLP_MSPE_CVA.band_lo.(f) = prctile(S2_81_hh_actual.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S2_81_hh_actual.aLP_MSPE_CVA.band_hi.(f) = prctile(S2_81_hh_actual.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S2_81_hh_actual.aLP_RS_GMA_girf.irf.(f) = S2_81_hh_actual.aLP_RS.alpha.(f) .* S2_81_hh_actual.GMA_lps.irf.(f)...
    + (1-S2_81_hh_actual.aLP_RS.alpha.(f)) .* S2_81_hh_actual.GMA_vars_girf.irf.(f); 
    S2_81_hh_actual.aLP_RS_GMA_girf.mu.(f) = mean(S2_81_hh_actual.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.aLP_RS_GMA_girf.band_lo.(f) = prctile(S2_81_hh_actual.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S2_81_hh_actual.aLP_RS_GMA_girf.band_hi.(f) = prctile(S2_81_hh_actual.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_81_hh_actual.aLP_RS_CVA_girf.irf.(f) = S2_81_hh_actual.aLP_RS.alpha.(f) .* S2_81_hh_actual.CVA_lps.irf.(f)...
    + (1-S2_81_hh_actual.aLP_RS.alpha.(f)) .* S2_81_hh_actual.CVA_vars_girf.irf.(f); 
    S2_81_hh_actual.aLP_RS_CVA_girf.mu.(f) = mean(S2_81_hh_actual.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.aLP_RS_CVA_girf.band_lo.(f) = prctile(S2_81_hh_actual.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S2_81_hh_actual.aLP_RS_CVA_girf.band_hi.(f) = prctile(S2_81_hh_actual.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S2_81_hh_actual.aLP_MSPE_GMA_girf.irf.(f) = S2_81_hh_actual.aLP_MSPE.alpha.(f) .* S2_81_hh_actual.GMA_lps.irf.(f)...
    + (1-S2_81_hh_actual.aLP_MSPE.alpha.(f)) .* S2_81_hh_actual.GMA_vars_girf.irf.(f); 
    S2_81_hh_actual.aLP_MSPE_GMA_girf.mu.(f) = mean(S2_81_hh_actual.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S2_81_hh_actual.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S2_81_hh_actual.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S2_81_hh_actual.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_81_hh_actual.aLP_MSPE_CVA_girf.irf.(f) = S2_81_hh_actual.aLP_MSPE.alpha.(f) .* S2_81_hh_actual.CVA_lps.irf.(f)...
    + (1-S2_81_hh_actual.aLP_MSPE.alpha.(f)) .* S2_81_hh_actual.CVA_vars_girf.irf.(f); 
    S2_81_hh_actual.aLP_MSPE_CVA_girf.mu.(f) = mean(S2_81_hh_actual.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_actual.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S2_81_hh_actual.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S2_81_hh_actual.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S2_81_hh_actual.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S2_81_hh_actual','-append')
%% ------------------------------------------------------
% 2-3. SPF IQR from 1981
% S2_81_pro_actual
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\actual\lag12\Bloom_S2_81_pro_actual.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S2_81_pro_actual.(names{i}).irf.(f) = irfpaths.(names{i}).S2_81_pro.(f);
        S2_81_pro_actual.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S2_81_pro.(f),2);
        S2_81_pro_actual.(names{i}).band_hi.(f) = band_hi.(names{i}).S2_81_pro.(f);
        S2_81_pro_actual.(names{i}).band_lo.(f) = band_lo.(names{i}).S2_81_pro.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S2_81_pro.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S2_81_pro_actual.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_81_pro.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S2_81_pro_actual.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_81_pro.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S2_81_pro_actual.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S2_81_pro.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S2_81_pro_actual.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S2_81_pro.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S2_81_pro_actual.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S2_81_pro.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S2_81_pro_actual.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_actual.GMA_lps.irf.(f) = wirf_val;% H x B
    S2_81_pro_actual.GMA_lps.mu.(f) = mean(S2_81_pro_actual.GMA_lps.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_actual.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S2_81_pro_actual.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_actual.CVA_lps.irf.(f) = wirf_val;% H x B
    S2_81_pro_actual.CVA_lps.mu.(f) = mean(S2_81_pro_actual.CVA_lps.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_actual.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S2_81_pro.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_81_pro_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_actual.GMA_vars.irf.(f) = wirf_val;% H x B
    S2_81_pro_actual.GMA_vars.mu.(f) = mean(S2_81_pro_actual.GMA_vars.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_actual.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_81_pro_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_actual.CVA_vars.irf.(f) = wirf_val;
    S2_81_pro_actual.CVA_vars.mu.(f) = mean(S2_81_pro_actual.CVA_vars.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_actual.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S2_81_pro.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_81_pro_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_actual.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S2_81_pro_actual.GMA_vars_girf.mu.(f) = mean(S2_81_pro_actual.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_actual.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_81_pro_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_actual.CVA_vars_girf.irf.(f) = wirf_val;
    S2_81_pro_actual.CVA_vars_girf.mu.(f) = mean(S2_81_pro_actual.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_actual.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S2_81_pro.state0_vol_empl = tss.S2_81_pro.state0_lps_empl;
    tss.S2_81_pro.state0_ff_empl = tss.S2_81_pro.state0_lps_empl;
    tss.S2_81_pro.state1_vol_empl = tss.S2_81_pro.state1_lps_empl;
    tss.S2_81_pro.state1_ff_empl = tss.S2_81_pro.state1_lps_empl;
    tss.S2_81_pro.state0_vol_ip = tss.S2_81_pro.state0_lps_ip;
    tss.S2_81_pro.state0_ff_ip = tss.S2_81_pro.state0_lps_ip;
    tss.S2_81_pro.state1_vol_ip = tss.S2_81_pro.state1_lps_ip;
    tss.S2_81_pro.state1_ff_ip = tss.S2_81_pro.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S2_81_pro.(f) ./ tss.S2_81_pro.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S2_81_pro.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S2_81_pro.state0_vol_empl = tss.S2_81_pro.state0_vars_empl;
    tss.S2_81_pro.state0_ff_empl = tss.S2_81_pro.state0_vars_empl;
    tss.S2_81_pro.state1_vol_empl = tss.S2_81_pro.state1_vars_empl;
    tss.S2_81_pro.state1_ff_empl = tss.S2_81_pro.state1_vars_empl;
    tss.S2_81_pro.state0_vol_ip = tss.S2_81_pro.state0_vars_ip;
    tss.S2_81_pro.state0_ff_ip = tss.S2_81_pro.state0_vars_ip;
    tss.S2_81_pro.state1_vol_ip = tss.S2_81_pro.state1_vars_ip;
    tss.S2_81_pro.state1_ff_ip = tss.S2_81_pro.state1_vars_ip;
    tss.S2_81_pro.linear_vol_empl = tss.S2_81_pro.linear_vars_empl;
    tss.S2_81_pro.linear_vol_ip = tss.S2_81_pro.linear_vars_ip;
    tss.S2_81_pro.linear_ff_empl = tss.S2_81_pro.linear_vars_empl;
    tss.S2_81_pro.linear_ff_ip = tss.S2_81_pro.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S2_81_pro.(f) ./ tss.S2_81_pro.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S2_81_pro.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S2_81_pro_actual.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S2_81_pro_actual.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S2_81_pro_actual.aLP_RS.alpha.(f) .* S2_81_pro_actual.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S2_81_pro_actual.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_81_pro_actual.GMA_vars.weight.(f)'; % H x 2
    S2_81_pro_actual.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S2_81_pro_actual.aLP_RS.alpha.(f) .* S2_81_pro_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_81_pro_actual.CVA_vars.weight.(f)';
    S2_81_pro_actual.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S2_81_pro_actual.aLP_MSPE.alpha.(f) .* S2_81_pro_actual.GMA_lps.weight.(f)';
    alpha_var = 1 - S2_81_pro_actual.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_81_pro_actual.GMA_vars.weight.(f)';
    S2_81_pro_actual.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S2_81_pro_actual.aLP_MSPE.alpha.(f) .* S2_81_pro_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_81_pro_actual.CVA_vars.weight.(f)';
    S2_81_pro_actual.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S2_81_pro_actual.aLP_RS_GMA.irf.(f) = S2_81_pro_actual.aLP_RS.alpha.(f) .* S2_81_pro_actual.GMA_lps.irf.(f)...
    + (1-S2_81_pro_actual.aLP_RS.alpha.(f)) .* S2_81_pro_actual.GMA_vars.irf.(f); 
    S2_81_pro_actual.aLP_RS_GMA.mu.(f) = mean(S2_81_pro_actual.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.aLP_RS_GMA.band_lo.(f) = prctile(S2_81_pro_actual.aLP_RS_GMA.irf.(f),pct_lo,2);
    S2_81_pro_actual.aLP_RS_GMA.band_hi.(f) = prctile(S2_81_pro_actual.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S2_81_pro_actual.aLP_RS_CVA.irf.(f) = S2_81_pro_actual.aLP_RS.alpha.(f) .* S2_81_pro_actual.CVA_lps.irf.(f)...
    + (1-S2_81_pro_actual.aLP_RS.alpha.(f)) .* S2_81_pro_actual.CVA_vars.irf.(f); 
    S2_81_pro_actual.aLP_RS_CVA.mu.(f) = mean(S2_81_pro_actual.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.aLP_RS_CVA.band_lo.(f) = prctile(S2_81_pro_actual.aLP_RS_CVA.irf.(f),pct_lo,2);
    S2_81_pro_actual.aLP_RS_CVA.band_hi.(f) = prctile(S2_81_pro_actual.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S2_81_pro_actual.aLP_MSPE_GMA.irf.(f) = S2_81_pro_actual.aLP_MSPE.alpha.(f) .* S2_81_pro_actual.GMA_lps.irf.(f)...
    + (1-S2_81_pro_actual.aLP_MSPE.alpha.(f)) .* S2_81_pro_actual.GMA_vars.irf.(f); 
    S2_81_pro_actual.aLP_MSPE_GMA.mu.(f) = mean(S2_81_pro_actual.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.aLP_MSPE_GMA.band_lo.(f) = prctile(S2_81_pro_actual.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S2_81_pro_actual.aLP_MSPE_GMA.band_hi.(f) = prctile(S2_81_pro_actual.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S2_81_pro_actual.aLP_MSPE_CVA.irf.(f) = S2_81_pro_actual.aLP_MSPE.alpha.(f) .* S2_81_pro_actual.CVA_lps.irf.(f)...
    + (1-S2_81_pro_actual.aLP_MSPE.alpha.(f)) .* S2_81_pro_actual.CVA_vars.irf.(f); 
    S2_81_pro_actual.aLP_MSPE_CVA.mu.(f) = mean(S2_81_pro_actual.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.aLP_MSPE_CVA.band_lo.(f) = prctile(S2_81_pro_actual.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S2_81_pro_actual.aLP_MSPE_CVA.band_hi.(f) = prctile(S2_81_pro_actual.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S2_81_pro_actual.aLP_RS_GMA_girf.irf.(f) = S2_81_pro_actual.aLP_RS.alpha.(f) .* S2_81_pro_actual.GMA_lps.irf.(f)...
    + (1-S2_81_pro_actual.aLP_RS.alpha.(f)) .* S2_81_pro_actual.GMA_vars_girf.irf.(f); 
    S2_81_pro_actual.aLP_RS_GMA_girf.mu.(f) = mean(S2_81_pro_actual.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.aLP_RS_GMA_girf.band_lo.(f) = prctile(S2_81_pro_actual.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S2_81_pro_actual.aLP_RS_GMA_girf.band_hi.(f) = prctile(S2_81_pro_actual.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_81_pro_actual.aLP_RS_CVA_girf.irf.(f) = S2_81_pro_actual.aLP_RS.alpha.(f) .* S2_81_pro_actual.CVA_lps.irf.(f)...
    + (1-S2_81_pro_actual.aLP_RS.alpha.(f)) .* S2_81_pro_actual.CVA_vars_girf.irf.(f); 
    S2_81_pro_actual.aLP_RS_CVA_girf.mu.(f) = mean(S2_81_pro_actual.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.aLP_RS_CVA_girf.band_lo.(f) = prctile(S2_81_pro_actual.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S2_81_pro_actual.aLP_RS_CVA_girf.band_hi.(f) = prctile(S2_81_pro_actual.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S2_81_pro_actual.aLP_MSPE_GMA_girf.irf.(f) = S2_81_pro_actual.aLP_MSPE.alpha.(f) .* S2_81_pro_actual.GMA_lps.irf.(f)...
    + (1-S2_81_pro_actual.aLP_MSPE.alpha.(f)) .* S2_81_pro_actual.GMA_vars_girf.irf.(f); 
    S2_81_pro_actual.aLP_MSPE_GMA_girf.mu.(f) = mean(S2_81_pro_actual.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S2_81_pro_actual.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S2_81_pro_actual.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S2_81_pro_actual.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_81_pro_actual.aLP_MSPE_CVA_girf.irf.(f) = S2_81_pro_actual.aLP_MSPE.alpha.(f) .* S2_81_pro_actual.CVA_lps.irf.(f)...
    + (1-S2_81_pro_actual.aLP_MSPE.alpha.(f)) .* S2_81_pro_actual.CVA_vars_girf.irf.(f); 
    S2_81_pro_actual.aLP_MSPE_CVA_girf.mu.(f) = mean(S2_81_pro_actual.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_actual.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S2_81_pro_actual.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S2_81_pro_actual.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S2_81_pro_actual.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S2_81_pro_actual','-append')
%% -------------------------------------------------------
% 3. Low vs. High inflation volatility
% 3-1. Fed target vs actual inflation rate by abolute
% S3_abs_actual
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\actual\lag12\Bloom_S3_abs_actual.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S3_abs_actual.(names{i}).irf.(f) = irfpaths.(names{i}).S3_abs.(f);
        S3_abs_actual.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S3_abs.(f),2);
        S3_abs_actual.(names{i}).band_hi.(f) = band_hi.(names{i}).S3_abs.(f);
        S3_abs_actual.(names{i}).band_lo.(f) = band_lo.(names{i}).S3_abs.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S3_abs.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S3_abs_actual.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_abs.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S3_abs_actual.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_abs.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S3_abs_actual.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S3_abs.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S3_abs_actual.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S3_abs.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S3_abs_actual.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S3_abs.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S3_abs_actual.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_actual.GMA_lps.irf.(f) = wirf_val;% H x B
    S3_abs_actual.GMA_lps.mu.(f) = mean(S3_abs_actual.GMA_lps.irf.(f), 2, "omitmissing");
    S3_abs_actual.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_actual.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S3_abs_actual.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_actual.CVA_lps.irf.(f) = wirf_val;% H x B
    S3_abs_actual.CVA_lps.mu.(f) = mean(S3_abs_actual.CVA_lps.irf.(f), 2, "omitmissing");
    S3_abs_actual.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_actual.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S3_abs.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_abs_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_actual.GMA_vars.irf.(f) = wirf_val;% H x B
    S3_abs_actual.GMA_vars.mu.(f) = mean(S3_abs_actual.GMA_vars.irf.(f), 2, "omitmissing");
    S3_abs_actual.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_actual.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_abs_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_actual.CVA_vars.irf.(f) = wirf_val;
    S3_abs_actual.CVA_vars.mu.(f) = mean(S3_abs_actual.CVA_vars.irf.(f), 2, "omitmissing");
    S3_abs_actual.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_actual.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S3_abs.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_abs_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_actual.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S3_abs_actual.GMA_vars_girf.mu.(f) = mean(S3_abs_actual.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S3_abs_actual.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_actual.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_abs_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_actual.CVA_vars_girf.irf.(f) = wirf_val;
    S3_abs_actual.CVA_vars_girf.mu.(f) = mean(S3_abs_actual.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S3_abs_actual.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_actual.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S3_abs.state0_vol_empl = tss.S3_abs.state0_lps_empl;
    tss.S3_abs.state0_ff_empl = tss.S3_abs.state0_lps_empl;
    tss.S3_abs.state1_vol_empl = tss.S3_abs.state1_lps_empl;
    tss.S3_abs.state1_ff_empl = tss.S3_abs.state1_lps_empl;
    tss.S3_abs.state0_vol_ip = tss.S3_abs.state0_lps_ip;
    tss.S3_abs.state0_ff_ip = tss.S3_abs.state0_lps_ip;
    tss.S3_abs.state1_vol_ip = tss.S3_abs.state1_lps_ip;
    tss.S3_abs.state1_ff_ip = tss.S3_abs.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S3_abs.(f) ./ tss.S3_abs.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S3_abs.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S3_abs.state0_vol_empl = tss.S3_abs.state0_vars_empl;
    tss.S3_abs.state0_ff_empl = tss.S3_abs.state0_vars_empl;
    tss.S3_abs.state1_vol_empl = tss.S3_abs.state1_vars_empl;
    tss.S3_abs.state1_ff_empl = tss.S3_abs.state1_vars_empl;
    tss.S3_abs.state0_vol_ip = tss.S3_abs.state0_vars_ip;
    tss.S3_abs.state0_ff_ip = tss.S3_abs.state0_vars_ip;
    tss.S3_abs.state1_vol_ip = tss.S3_abs.state1_vars_ip;
    tss.S3_abs.state1_ff_ip = tss.S3_abs.state1_vars_ip;
    tss.S3_abs.linear_vol_empl = tss.S3_abs.linear_vars_empl;
    tss.S3_abs.linear_vol_ip = tss.S3_abs.linear_vars_ip;
    tss.S3_abs.linear_ff_empl = tss.S3_abs.linear_vars_empl;
    tss.S3_abs.linear_ff_ip = tss.S3_abs.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S3_abs.(f) ./ tss.S3_abs.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S3_abs.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S3_abs_actual.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S3_abs_actual.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S3_abs_actual.aLP_RS.alpha.(f) .* S3_abs_actual.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S3_abs_actual.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_abs_actual.GMA_vars.weight.(f)'; % H x 2
    S3_abs_actual.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S3_abs_actual.aLP_RS.alpha.(f) .* S3_abs_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_abs_actual.CVA_vars.weight.(f)';
    S3_abs_actual.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S3_abs_actual.aLP_MSPE.alpha.(f) .* S3_abs_actual.GMA_lps.weight.(f)';
    alpha_var = 1 - S3_abs_actual.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_abs_actual.GMA_vars.weight.(f)';
    S3_abs_actual.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S3_abs_actual.aLP_MSPE.alpha.(f) .* S3_abs_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_abs_actual.CVA_vars.weight.(f)';
    S3_abs_actual.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S3_abs_actual.aLP_RS_GMA.irf.(f) = S3_abs_actual.aLP_RS.alpha.(f) .* S3_abs_actual.GMA_lps.irf.(f)...
    + (1-S3_abs_actual.aLP_RS.alpha.(f)) .* S3_abs_actual.GMA_vars.irf.(f); 
    S3_abs_actual.aLP_RS_GMA.mu.(f) = mean(S3_abs_actual.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S3_abs_actual.aLP_RS_GMA.band_lo.(f) = prctile(S3_abs_actual.aLP_RS_GMA.irf.(f),pct_lo,2);
    S3_abs_actual.aLP_RS_GMA.band_hi.(f) = prctile(S3_abs_actual.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S3_abs_actual.aLP_RS_CVA.irf.(f) = S3_abs_actual.aLP_RS.alpha.(f) .* S3_abs_actual.CVA_lps.irf.(f)...
    + (1-S3_abs_actual.aLP_RS.alpha.(f)) .* S3_abs_actual.CVA_vars.irf.(f); 
    S3_abs_actual.aLP_RS_CVA.mu.(f) = mean(S3_abs_actual.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S3_abs_actual.aLP_RS_CVA.band_lo.(f) = prctile(S3_abs_actual.aLP_RS_CVA.irf.(f),pct_lo,2);
    S3_abs_actual.aLP_RS_CVA.band_hi.(f) = prctile(S3_abs_actual.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S3_abs_actual.aLP_MSPE_GMA.irf.(f) = S3_abs_actual.aLP_MSPE.alpha.(f) .* S3_abs_actual.GMA_lps.irf.(f)...
    + (1-S3_abs_actual.aLP_MSPE.alpha.(f)) .* S3_abs_actual.GMA_vars.irf.(f); 
    S3_abs_actual.aLP_MSPE_GMA.mu.(f) = mean(S3_abs_actual.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S3_abs_actual.aLP_MSPE_GMA.band_lo.(f) = prctile(S3_abs_actual.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S3_abs_actual.aLP_MSPE_GMA.band_hi.(f) = prctile(S3_abs_actual.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S3_abs_actual.aLP_MSPE_CVA.irf.(f) = S3_abs_actual.aLP_MSPE.alpha.(f) .* S3_abs_actual.CVA_lps.irf.(f)...
    + (1-S3_abs_actual.aLP_MSPE.alpha.(f)) .* S3_abs_actual.CVA_vars.irf.(f); 
    S3_abs_actual.aLP_MSPE_CVA.mu.(f) = mean(S3_abs_actual.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S3_abs_actual.aLP_MSPE_CVA.band_lo.(f) = prctile(S3_abs_actual.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S3_abs_actual.aLP_MSPE_CVA.band_hi.(f) = prctile(S3_abs_actual.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S3_abs_actual.aLP_RS_GMA_girf.irf.(f) = S3_abs_actual.aLP_RS.alpha.(f) .* S3_abs_actual.GMA_lps.irf.(f)...
    + (1-S3_abs_actual.aLP_RS.alpha.(f)) .* S3_abs_actual.GMA_vars_girf.irf.(f); 
    S3_abs_actual.aLP_RS_GMA_girf.mu.(f) = mean(S3_abs_actual.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S3_abs_actual.aLP_RS_GMA_girf.band_lo.(f) = prctile(S3_abs_actual.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S3_abs_actual.aLP_RS_GMA_girf.band_hi.(f) = prctile(S3_abs_actual.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_abs_actual.aLP_RS_CVA_girf.irf.(f) = S3_abs_actual.aLP_RS.alpha.(f) .* S3_abs_actual.CVA_lps.irf.(f)...
    + (1-S3_abs_actual.aLP_RS.alpha.(f)) .* S3_abs_actual.CVA_vars_girf.irf.(f); 
    S3_abs_actual.aLP_RS_CVA_girf.mu.(f) = mean(S3_abs_actual.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S3_abs_actual.aLP_RS_CVA_girf.band_lo.(f) = prctile(S3_abs_actual.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S3_abs_actual.aLP_RS_CVA_girf.band_hi.(f) = prctile(S3_abs_actual.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S3_abs_actual.aLP_MSPE_GMA_girf.irf.(f) = S3_abs_actual.aLP_MSPE.alpha.(f) .* S3_abs_actual.GMA_lps.irf.(f)...
    + (1-S3_abs_actual.aLP_MSPE.alpha.(f)) .* S3_abs_actual.GMA_vars_girf.irf.(f); 
    S3_abs_actual.aLP_MSPE_GMA_girf.mu.(f) = mean(S3_abs_actual.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S3_abs_actual.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S3_abs_actual.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S3_abs_actual.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S3_abs_actual.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_abs_actual.aLP_MSPE_CVA_girf.irf.(f) = S3_abs_actual.aLP_MSPE.alpha.(f) .* S3_abs_actual.CVA_lps.irf.(f)...
    + (1-S3_abs_actual.aLP_MSPE.alpha.(f)) .* S3_abs_actual.CVA_vars_girf.irf.(f); 
    S3_abs_actual.aLP_MSPE_CVA_girf.mu.(f) = mean(S3_abs_actual.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S3_abs_actual.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S3_abs_actual.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S3_abs_actual.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S3_abs_actual.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S3_abs_actual','-append')
%% ------------------------------------------------------
% 3-2. Fed target vs actual inflation rate by 75th percentile
% S3_rel_actual
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\actual\lag12\Bloom_S3_rel_actual.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S3_rel_actual.(names{i}).irf.(f) = irfpaths.(names{i}).S3_rel.(f);
        S3_rel_actual.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S3_rel.(f),2);
        S3_rel_actual.(names{i}).band_hi.(f) = band_hi.(names{i}).S3_rel.(f);
        S3_rel_actual.(names{i}).band_lo.(f) = band_lo.(names{i}).S3_rel.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S3_rel.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S3_rel_actual.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_rel.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S3_rel_actual.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_rel.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S3_rel_actual.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S3_rel.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S3_rel_actual.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S3_rel.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S3_rel_actual.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S3_rel.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S3_rel_actual.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_actual.GMA_lps.irf.(f) = wirf_val;% H x B
    S3_rel_actual.GMA_lps.mu.(f) = mean(S3_rel_actual.GMA_lps.irf.(f), 2, "omitmissing");
    S3_rel_actual.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_actual.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S3_rel_actual.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_actual.CVA_lps.irf.(f) = wirf_val;% H x B
    S3_rel_actual.CVA_lps.mu.(f) = mean(S3_rel_actual.CVA_lps.irf.(f), 2, "omitmissing");
    S3_rel_actual.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_actual.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S3_rel.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_rel_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_actual.GMA_vars.irf.(f) = wirf_val;% H x B
    S3_rel_actual.GMA_vars.mu.(f) = mean(S3_rel_actual.GMA_vars.irf.(f), 2, "omitmissing");
    S3_rel_actual.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_actual.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_rel_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_actual.CVA_vars.irf.(f) = wirf_val;
    S3_rel_actual.CVA_vars.mu.(f) = mean(S3_rel_actual.CVA_vars.irf.(f), 2, "omitmissing");
    S3_rel_actual.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_actual.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S3_rel.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_rel_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_actual.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S3_rel_actual.GMA_vars_girf.mu.(f) = mean(S3_rel_actual.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S3_rel_actual.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_actual.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_rel_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_actual.CVA_vars_girf.irf.(f) = wirf_val;
    S3_rel_actual.CVA_vars_girf.mu.(f) = mean(S3_rel_actual.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S3_rel_actual.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_actual.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S3_rel.state0_vol_empl = tss.S3_rel.state0_lps_empl;
    tss.S3_rel.state0_ff_empl = tss.S3_rel.state0_lps_empl;
    tss.S3_rel.state1_vol_empl = tss.S3_rel.state1_lps_empl;
    tss.S3_rel.state1_ff_empl = tss.S3_rel.state1_lps_empl;
    tss.S3_rel.state0_vol_ip = tss.S3_rel.state0_lps_ip;
    tss.S3_rel.state0_ff_ip = tss.S3_rel.state0_lps_ip;
    tss.S3_rel.state1_vol_ip = tss.S3_rel.state1_lps_ip;
    tss.S3_rel.state1_ff_ip = tss.S3_rel.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S3_rel.(f) ./ tss.S3_rel.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S3_rel.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S3_rel.state0_vol_empl = tss.S3_rel.state0_vars_empl;
    tss.S3_rel.state0_ff_empl = tss.S3_rel.state0_vars_empl;
    tss.S3_rel.state1_vol_empl = tss.S3_rel.state1_vars_empl;
    tss.S3_rel.state1_ff_empl = tss.S3_rel.state1_vars_empl;
    tss.S3_rel.state0_vol_ip = tss.S3_rel.state0_vars_ip;
    tss.S3_rel.state0_ff_ip = tss.S3_rel.state0_vars_ip;
    tss.S3_rel.state1_vol_ip = tss.S3_rel.state1_vars_ip;
    tss.S3_rel.state1_ff_ip = tss.S3_rel.state1_vars_ip;
    tss.S3_rel.linear_vol_empl = tss.S3_rel.linear_vars_empl;
    tss.S3_rel.linear_vol_ip = tss.S3_rel.linear_vars_ip;
    tss.S3_rel.linear_ff_empl = tss.S3_rel.linear_vars_empl;
    tss.S3_rel.linear_ff_ip = tss.S3_rel.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S3_rel.(f) ./ tss.S3_rel.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S3_rel.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S3_rel_actual.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S3_rel_actual.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S3_rel_actual.aLP_RS.alpha.(f) .* S3_rel_actual.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S3_rel_actual.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_rel_actual.GMA_vars.weight.(f)'; % H x 2
    S3_rel_actual.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S3_rel_actual.aLP_RS.alpha.(f) .* S3_rel_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_rel_actual.CVA_vars.weight.(f)';
    S3_rel_actual.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S3_rel_actual.aLP_MSPE.alpha.(f) .* S3_rel_actual.GMA_lps.weight.(f)';
    alpha_var = 1 - S3_rel_actual.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_rel_actual.GMA_vars.weight.(f)';
    S3_rel_actual.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S3_rel_actual.aLP_MSPE.alpha.(f) .* S3_rel_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_rel_actual.CVA_vars.weight.(f)';
    S3_rel_actual.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S3_rel_actual.aLP_RS_GMA.irf.(f) = S3_rel_actual.aLP_RS.alpha.(f) .* S3_rel_actual.GMA_lps.irf.(f)...
    + (1-S3_rel_actual.aLP_RS.alpha.(f)) .* S3_rel_actual.GMA_vars.irf.(f); 
    S3_rel_actual.aLP_RS_GMA.mu.(f) = mean(S3_rel_actual.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S3_rel_actual.aLP_RS_GMA.band_lo.(f) = prctile(S3_rel_actual.aLP_RS_GMA.irf.(f),pct_lo,2);
    S3_rel_actual.aLP_RS_GMA.band_hi.(f) = prctile(S3_rel_actual.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S3_rel_actual.aLP_RS_CVA.irf.(f) = S3_rel_actual.aLP_RS.alpha.(f) .* S3_rel_actual.CVA_lps.irf.(f)...
    + (1-S3_rel_actual.aLP_RS.alpha.(f)) .* S3_rel_actual.CVA_vars.irf.(f); 
    S3_rel_actual.aLP_RS_CVA.mu.(f) = mean(S3_rel_actual.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S3_rel_actual.aLP_RS_CVA.band_lo.(f) = prctile(S3_rel_actual.aLP_RS_CVA.irf.(f),pct_lo,2);
    S3_rel_actual.aLP_RS_CVA.band_hi.(f) = prctile(S3_rel_actual.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S3_rel_actual.aLP_MSPE_GMA.irf.(f) = S3_rel_actual.aLP_MSPE.alpha.(f) .* S3_rel_actual.GMA_lps.irf.(f)...
    + (1-S3_rel_actual.aLP_MSPE.alpha.(f)) .* S3_rel_actual.GMA_vars.irf.(f); 
    S3_rel_actual.aLP_MSPE_GMA.mu.(f) = mean(S3_rel_actual.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S3_rel_actual.aLP_MSPE_GMA.band_lo.(f) = prctile(S3_rel_actual.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S3_rel_actual.aLP_MSPE_GMA.band_hi.(f) = prctile(S3_rel_actual.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S3_rel_actual.aLP_MSPE_CVA.irf.(f) = S3_rel_actual.aLP_MSPE.alpha.(f) .* S3_rel_actual.CVA_lps.irf.(f)...
    + (1-S3_rel_actual.aLP_MSPE.alpha.(f)) .* S3_rel_actual.CVA_vars.irf.(f); 
    S3_rel_actual.aLP_MSPE_CVA.mu.(f) = mean(S3_rel_actual.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S3_rel_actual.aLP_MSPE_CVA.band_lo.(f) = prctile(S3_rel_actual.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S3_rel_actual.aLP_MSPE_CVA.band_hi.(f) = prctile(S3_rel_actual.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S3_rel_actual.aLP_RS_GMA_girf.irf.(f) = S3_rel_actual.aLP_RS.alpha.(f) .* S3_rel_actual.GMA_lps.irf.(f)...
    + (1-S3_rel_actual.aLP_RS.alpha.(f)) .* S3_rel_actual.GMA_vars_girf.irf.(f); 
    S3_rel_actual.aLP_RS_GMA_girf.mu.(f) = mean(S3_rel_actual.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S3_rel_actual.aLP_RS_GMA_girf.band_lo.(f) = prctile(S3_rel_actual.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S3_rel_actual.aLP_RS_GMA_girf.band_hi.(f) = prctile(S3_rel_actual.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_rel_actual.aLP_RS_CVA_girf.irf.(f) = S3_rel_actual.aLP_RS.alpha.(f) .* S3_rel_actual.CVA_lps.irf.(f)...
    + (1-S3_rel_actual.aLP_RS.alpha.(f)) .* S3_rel_actual.CVA_vars_girf.irf.(f); 
    S3_rel_actual.aLP_RS_CVA_girf.mu.(f) = mean(S3_rel_actual.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S3_rel_actual.aLP_RS_CVA_girf.band_lo.(f) = prctile(S3_rel_actual.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S3_rel_actual.aLP_RS_CVA_girf.band_hi.(f) = prctile(S3_rel_actual.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S3_rel_actual.aLP_MSPE_GMA_girf.irf.(f) = S3_rel_actual.aLP_MSPE.alpha.(f) .* S3_rel_actual.GMA_lps.irf.(f)...
    + (1-S3_rel_actual.aLP_MSPE.alpha.(f)) .* S3_rel_actual.GMA_vars_girf.irf.(f); 
    S3_rel_actual.aLP_MSPE_GMA_girf.mu.(f) = mean(S3_rel_actual.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S3_rel_actual.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S3_rel_actual.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S3_rel_actual.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S3_rel_actual.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_rel_actual.aLP_MSPE_CVA_girf.irf.(f) = S3_rel_actual.aLP_MSPE.alpha.(f) .* S3_rel_actual.CVA_lps.irf.(f)...
    + (1-S3_rel_actual.aLP_MSPE.alpha.(f)) .* S3_rel_actual.CVA_vars_girf.irf.(f); 
    S3_rel_actual.aLP_MSPE_CVA_girf.mu.(f) = mean(S3_rel_actual.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S3_rel_actual.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S3_rel_actual.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S3_rel_actual.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S3_rel_actual.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S3_rel_actual','-append')
%% ------------------------------------------------------
% 3-3. by rolling sigma 75th percentile
% S3_75thr_actual
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\actual\lag12\Bloom_S3_75thr_actual.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S3_75thr_actual.(names{i}).irf.(f) = irfpaths.(names{i}).S3_75thr.(f);
        S3_75thr_actual.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S3_75thr.(f),2);
        S3_75thr_actual.(names{i}).band_hi.(f) = band_hi.(names{i}).S3_75thr.(f);
        S3_75thr_actual.(names{i}).band_lo.(f) = band_lo.(names{i}).S3_75thr.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S3_75thr.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S3_75thr_actual.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_75thr.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S3_75thr_actual.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_75thr.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S3_75thr_actual.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S3_75thr.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S3_75thr_actual.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S3_75thr.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S3_75thr_actual.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S3_75thr.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S3_75thr_actual.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_actual.GMA_lps.irf.(f) = wirf_val;% H x B
    S3_75thr_actual.GMA_lps.mu.(f) = mean(S3_75thr_actual.GMA_lps.irf.(f), 2, "omitmissing");
    S3_75thr_actual.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_actual.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S3_75thr_actual.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_actual.CVA_lps.irf.(f) = wirf_val;% H x B
    S3_75thr_actual.CVA_lps.mu.(f) = mean(S3_75thr_actual.CVA_lps.irf.(f), 2, "omitmissing");
    S3_75thr_actual.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_actual.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S3_75thr.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_75thr_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_actual.GMA_vars.irf.(f) = wirf_val;% H x B
    S3_75thr_actual.GMA_vars.mu.(f) = mean(S3_75thr_actual.GMA_vars.irf.(f), 2, "omitmissing");
    S3_75thr_actual.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_actual.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_75thr_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_actual.CVA_vars.irf.(f) = wirf_val;
    S3_75thr_actual.CVA_vars.mu.(f) = mean(S3_75thr_actual.CVA_vars.irf.(f), 2, "omitmissing");
    S3_75thr_actual.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_actual.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S3_75thr.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_75thr_actual.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_actual.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S3_75thr_actual.GMA_vars_girf.mu.(f) = mean(S3_75thr_actual.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S3_75thr_actual.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_actual.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_75thr_actual.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_actual.CVA_vars_girf.irf.(f) = wirf_val;
    S3_75thr_actual.CVA_vars_girf.mu.(f) = mean(S3_75thr_actual.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S3_75thr_actual.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_actual.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S3_75thr.state0_vol_empl = tss.S3_75thr.state0_lps_empl;
    tss.S3_75thr.state0_ff_empl = tss.S3_75thr.state0_lps_empl;
    tss.S3_75thr.state1_vol_empl = tss.S3_75thr.state1_lps_empl;
    tss.S3_75thr.state1_ff_empl = tss.S3_75thr.state1_lps_empl;
    tss.S3_75thr.state0_vol_ip = tss.S3_75thr.state0_lps_ip;
    tss.S3_75thr.state0_ff_ip = tss.S3_75thr.state0_lps_ip;
    tss.S3_75thr.state1_vol_ip = tss.S3_75thr.state1_lps_ip;
    tss.S3_75thr.state1_ff_ip = tss.S3_75thr.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S3_75thr.(f) ./ tss.S3_75thr.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S3_75thr.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S3_75thr.state0_vol_empl = tss.S3_75thr.state0_vars_empl;
    tss.S3_75thr.state0_ff_empl = tss.S3_75thr.state0_vars_empl;
    tss.S3_75thr.state1_vol_empl = tss.S3_75thr.state1_vars_empl;
    tss.S3_75thr.state1_ff_empl = tss.S3_75thr.state1_vars_empl;
    tss.S3_75thr.state0_vol_ip = tss.S3_75thr.state0_vars_ip;
    tss.S3_75thr.state0_ff_ip = tss.S3_75thr.state0_vars_ip;
    tss.S3_75thr.state1_vol_ip = tss.S3_75thr.state1_vars_ip;
    tss.S3_75thr.state1_ff_ip = tss.S3_75thr.state1_vars_ip;
    tss.S3_75thr.linear_vol_empl = tss.S3_75thr.linear_vars_empl;
    tss.S3_75thr.linear_vol_ip = tss.S3_75thr.linear_vars_ip;
    tss.S3_75thr.linear_ff_empl = tss.S3_75thr.linear_vars_empl;
    tss.S3_75thr.linear_ff_ip = tss.S3_75thr.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S3_75thr.(f) ./ tss.S3_75thr.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S3_75thr.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S3_75thr_actual.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S3_75thr_actual.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S3_75thr_actual.aLP_RS.alpha.(f) .* S3_75thr_actual.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S3_75thr_actual.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_75thr_actual.GMA_vars.weight.(f)'; % H x 2
    S3_75thr_actual.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S3_75thr_actual.aLP_RS.alpha.(f) .* S3_75thr_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_75thr_actual.CVA_vars.weight.(f)';
    S3_75thr_actual.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S3_75thr_actual.aLP_MSPE.alpha.(f) .* S3_75thr_actual.GMA_lps.weight.(f)';
    alpha_var = 1 - S3_75thr_actual.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_75thr_actual.GMA_vars.weight.(f)';
    S3_75thr_actual.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S3_75thr_actual.aLP_MSPE.alpha.(f) .* S3_75thr_actual.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_75thr_actual.CVA_vars.weight.(f)';
    S3_75thr_actual.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S3_75thr_actual.aLP_RS_GMA.irf.(f) = S3_75thr_actual.aLP_RS.alpha.(f) .* S3_75thr_actual.GMA_lps.irf.(f)...
    + (1-S3_75thr_actual.aLP_RS.alpha.(f)) .* S3_75thr_actual.GMA_vars.irf.(f); 
    S3_75thr_actual.aLP_RS_GMA.mu.(f) = mean(S3_75thr_actual.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S3_75thr_actual.aLP_RS_GMA.band_lo.(f) = prctile(S3_75thr_actual.aLP_RS_GMA.irf.(f),pct_lo,2);
    S3_75thr_actual.aLP_RS_GMA.band_hi.(f) = prctile(S3_75thr_actual.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S3_75thr_actual.aLP_RS_CVA.irf.(f) = S3_75thr_actual.aLP_RS.alpha.(f) .* S3_75thr_actual.CVA_lps.irf.(f)...
    + (1-S3_75thr_actual.aLP_RS.alpha.(f)) .* S3_75thr_actual.CVA_vars.irf.(f); 
    S3_75thr_actual.aLP_RS_CVA.mu.(f) = mean(S3_75thr_actual.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S3_75thr_actual.aLP_RS_CVA.band_lo.(f) = prctile(S3_75thr_actual.aLP_RS_CVA.irf.(f),pct_lo,2);
    S3_75thr_actual.aLP_RS_CVA.band_hi.(f) = prctile(S3_75thr_actual.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S3_75thr_actual.aLP_MSPE_GMA.irf.(f) = S3_75thr_actual.aLP_MSPE.alpha.(f) .* S3_75thr_actual.GMA_lps.irf.(f)...
    + (1-S3_75thr_actual.aLP_MSPE.alpha.(f)) .* S3_75thr_actual.GMA_vars.irf.(f); 
    S3_75thr_actual.aLP_MSPE_GMA.mu.(f) = mean(S3_75thr_actual.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S3_75thr_actual.aLP_MSPE_GMA.band_lo.(f) = prctile(S3_75thr_actual.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S3_75thr_actual.aLP_MSPE_GMA.band_hi.(f) = prctile(S3_75thr_actual.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S3_75thr_actual.aLP_MSPE_CVA.irf.(f) = S3_75thr_actual.aLP_MSPE.alpha.(f) .* S3_75thr_actual.CVA_lps.irf.(f)...
    + (1-S3_75thr_actual.aLP_MSPE.alpha.(f)) .* S3_75thr_actual.CVA_vars.irf.(f); 
    S3_75thr_actual.aLP_MSPE_CVA.mu.(f) = mean(S3_75thr_actual.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S3_75thr_actual.aLP_MSPE_CVA.band_lo.(f) = prctile(S3_75thr_actual.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S3_75thr_actual.aLP_MSPE_CVA.band_hi.(f) = prctile(S3_75thr_actual.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S3_75thr_actual.aLP_RS_GMA_girf.irf.(f) = S3_75thr_actual.aLP_RS.alpha.(f) .* S3_75thr_actual.GMA_lps.irf.(f)...
    + (1-S3_75thr_actual.aLP_RS.alpha.(f)) .* S3_75thr_actual.GMA_vars_girf.irf.(f); 
    S3_75thr_actual.aLP_RS_GMA_girf.mu.(f) = mean(S3_75thr_actual.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S3_75thr_actual.aLP_RS_GMA_girf.band_lo.(f) = prctile(S3_75thr_actual.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S3_75thr_actual.aLP_RS_GMA_girf.band_hi.(f) = prctile(S3_75thr_actual.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_75thr_actual.aLP_RS_CVA_girf.irf.(f) = S3_75thr_actual.aLP_RS.alpha.(f) .* S3_75thr_actual.CVA_lps.irf.(f)...
    + (1-S3_75thr_actual.aLP_RS.alpha.(f)) .* S3_75thr_actual.CVA_vars_girf.irf.(f); 
    S3_75thr_actual.aLP_RS_CVA_girf.mu.(f) = mean(S3_75thr_actual.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S3_75thr_actual.aLP_RS_CVA_girf.band_lo.(f) = prctile(S3_75thr_actual.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S3_75thr_actual.aLP_RS_CVA_girf.band_hi.(f) = prctile(S3_75thr_actual.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S3_75thr_actual.aLP_MSPE_GMA_girf.irf.(f) = S3_75thr_actual.aLP_MSPE.alpha.(f) .* S3_75thr_actual.GMA_lps.irf.(f)...
    + (1-S3_75thr_actual.aLP_MSPE.alpha.(f)) .* S3_75thr_actual.GMA_vars_girf.irf.(f); 
    S3_75thr_actual.aLP_MSPE_GMA_girf.mu.(f) = mean(S3_75thr_actual.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S3_75thr_actual.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S3_75thr_actual.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S3_75thr_actual.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S3_75thr_actual.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_75thr_actual.aLP_MSPE_CVA_girf.irf.(f) = S3_75thr_actual.aLP_MSPE.alpha.(f) .* S3_75thr_actual.CVA_lps.irf.(f)...
    + (1-S3_75thr_actual.aLP_MSPE.alpha.(f)) .* S3_75thr_actual.CVA_vars_girf.irf.(f); 
    S3_75thr_actual.aLP_MSPE_CVA_girf.mu.(f) = mean(S3_75thr_actual.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S3_75thr_actual.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S3_75thr_actual.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S3_75thr_actual.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S3_75thr_actual.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S3_75thr_actual','-append')
%% ------------------------------------------------------
% 2nd shock = Indicator 
% -------------------------------------------------------
% 1. Good and bad inflation regimes 
% 1-1. the rolling window (60 months) correlation between inflation (yoy) and
% consumption growth (actual)
% S1_60_id_indcator
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\indicator\lag12\Bloom_S1_60_id_indicator.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S1_60_id_indicator.(names{i}).irf.(f) = irfpaths.(names{i}).S1_60_id.(f);
        S1_60_id_indicator.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S1_60_id.(f),2);
        S1_60_id_indicator.(names{i}).band_hi.(f) = band_hi.(names{i}).S1_60_id.(f);
        S1_60_id_indicator.(names{i}).band_lo.(f) = band_lo.(names{i}).S1_60_id.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S1_60_id.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S1_60_id_indicator.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_60_id.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S1_60_id_indicator.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_60_id.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S1_60_id_indicator.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S1_60_id.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S1_60_id_indicator.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S1_60_id.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S1_60_id_indicator.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S1_60_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S1_60_id_indicator.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_indicator.GMA_lps.irf.(f) = wirf_val;% H x B
    S1_60_id_indicator.GMA_lps.mu.(f) = mean(S1_60_id_indicator.GMA_lps.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_indicator.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S1_60_id_indicator.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_indicator.CVA_lps.irf.(f) = wirf_val;% H x B
    S1_60_id_indicator.CVA_lps.mu.(f) = mean(S1_60_id_indicator.CVA_lps.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_indicator.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S1_60_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_60_id_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_indicator.GMA_vars.irf.(f) = wirf_val;% H x B
    S1_60_id_indicator.GMA_vars.mu.(f) = mean(S1_60_id_indicator.GMA_vars.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_indicator.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_60_id_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_indicator.CVA_vars.irf.(f) = wirf_val;
    S1_60_id_indicator.CVA_vars.mu.(f) = mean(S1_60_id_indicator.CVA_vars.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_indicator.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S1_60_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_60_id_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_indicator.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S1_60_id_indicator.GMA_vars_girf.mu.(f) = mean(S1_60_id_indicator.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_indicator.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_60_id_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_id_indicator.CVA_vars_girf.irf.(f) = wirf_val;
    S1_60_id_indicator.CVA_vars_girf.mu.(f) = mean(S1_60_id_indicator.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_id_indicator.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S1_60_id.state0_vol_empl = tss.S1_60_id.state0_lps_empl;
    tss.S1_60_id.state0_ff_empl = tss.S1_60_id.state0_lps_empl;
    tss.S1_60_id.state1_vol_empl = tss.S1_60_id.state1_lps_empl;
    tss.S1_60_id.state1_ff_empl = tss.S1_60_id.state1_lps_empl;
    tss.S1_60_id.state0_vol_ip = tss.S1_60_id.state0_lps_ip;
    tss.S1_60_id.state0_ff_ip = tss.S1_60_id.state0_lps_ip;
    tss.S1_60_id.state1_vol_ip = tss.S1_60_id.state1_lps_ip;
    tss.S1_60_id.state1_ff_ip = tss.S1_60_id.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S1_60_id.(f) ./ tss.S1_60_id.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S1_60_id.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S1_60_id.state0_vol_empl = tss.S1_60_id.state0_vars_empl;
    tss.S1_60_id.state0_ff_empl = tss.S1_60_id.state0_vars_empl;
    tss.S1_60_id.state1_vol_empl = tss.S1_60_id.state1_vars_empl;
    tss.S1_60_id.state1_ff_empl = tss.S1_60_id.state1_vars_empl;
    tss.S1_60_id.state0_vol_ip = tss.S1_60_id.state0_vars_ip;
    tss.S1_60_id.state0_ff_ip = tss.S1_60_id.state0_vars_ip;
    tss.S1_60_id.state1_vol_ip = tss.S1_60_id.state1_vars_ip;
    tss.S1_60_id.state1_ff_ip = tss.S1_60_id.state1_vars_ip;
    tss.S1_60_id.linear_vol_empl = tss.S1_60_id.linear_vars_empl;
    tss.S1_60_id.linear_vol_ip = tss.S1_60_id.linear_vars_ip;
    tss.S1_60_id.linear_ff_empl = tss.S1_60_id.linear_vars_empl;
    tss.S1_60_id.linear_ff_ip = tss.S1_60_id.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S1_60_id.(f) ./ tss.S1_60_id.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S1_60_id.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S1_60_id_indicator.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S1_60_id_indicator.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S1_60_id_indicator.aLP_RS.alpha.(f) .* S1_60_id_indicator.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S1_60_id_indicator.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_60_id_indicator.GMA_vars.weight.(f)'; % H x 2
    S1_60_id_indicator.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S1_60_id_indicator.aLP_RS.alpha.(f) .* S1_60_id_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_60_id_indicator.CVA_vars.weight.(f)';
    S1_60_id_indicator.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S1_60_id_indicator.aLP_MSPE.alpha.(f) .* S1_60_id_indicator.GMA_lps.weight.(f)';
    alpha_var = 1 - S1_60_id_indicator.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_60_id_indicator.GMA_vars.weight.(f)';
    S1_60_id_indicator.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S1_60_id_indicator.aLP_MSPE.alpha.(f) .* S1_60_id_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_60_id_indicator.CVA_vars.weight.(f)';
    S1_60_id_indicator.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S1_60_id_indicator.aLP_RS_GMA.irf.(f) = S1_60_id_indicator.aLP_RS.alpha.(f) .* S1_60_id_indicator.GMA_lps.irf.(f)...
    + (1-S1_60_id_indicator.aLP_RS.alpha.(f)) .* S1_60_id_indicator.GMA_vars.irf.(f); 
    S1_60_id_indicator.aLP_RS_GMA.mu.(f) = mean(S1_60_id_indicator.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.aLP_RS_GMA.band_lo.(f) = prctile(S1_60_id_indicator.aLP_RS_GMA.irf.(f),pct_lo,2);
    S1_60_id_indicator.aLP_RS_GMA.band_hi.(f) = prctile(S1_60_id_indicator.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S1_60_id_indicator.aLP_RS_CVA.irf.(f) = S1_60_id_indicator.aLP_RS.alpha.(f) .* S1_60_id_indicator.CVA_lps.irf.(f)...
    + (1-S1_60_id_indicator.aLP_RS.alpha.(f)) .* S1_60_id_indicator.CVA_vars.irf.(f); 
    S1_60_id_indicator.aLP_RS_CVA.mu.(f) = mean(S1_60_id_indicator.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.aLP_RS_CVA.band_lo.(f) = prctile(S1_60_id_indicator.aLP_RS_CVA.irf.(f),pct_lo,2);
    S1_60_id_indicator.aLP_RS_CVA.band_hi.(f) = prctile(S1_60_id_indicator.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S1_60_id_indicator.aLP_MSPE_GMA.irf.(f) = S1_60_id_indicator.aLP_MSPE.alpha.(f) .* S1_60_id_indicator.GMA_lps.irf.(f)...
    + (1-S1_60_id_indicator.aLP_MSPE.alpha.(f)) .* S1_60_id_indicator.GMA_vars.irf.(f); 
    S1_60_id_indicator.aLP_MSPE_GMA.mu.(f) = mean(S1_60_id_indicator.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.aLP_MSPE_GMA.band_lo.(f) = prctile(S1_60_id_indicator.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S1_60_id_indicator.aLP_MSPE_GMA.band_hi.(f) = prctile(S1_60_id_indicator.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S1_60_id_indicator.aLP_MSPE_CVA.irf.(f) = S1_60_id_indicator.aLP_MSPE.alpha.(f) .* S1_60_id_indicator.CVA_lps.irf.(f)...
    + (1-S1_60_id_indicator.aLP_MSPE.alpha.(f)) .* S1_60_id_indicator.CVA_vars.irf.(f); 
    S1_60_id_indicator.aLP_MSPE_CVA.mu.(f) = mean(S1_60_id_indicator.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.aLP_MSPE_CVA.band_lo.(f) = prctile(S1_60_id_indicator.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S1_60_id_indicator.aLP_MSPE_CVA.band_hi.(f) = prctile(S1_60_id_indicator.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S1_60_id_indicator.aLP_RS_GMA_girf.irf.(f) = S1_60_id_indicator.aLP_RS.alpha.(f) .* S1_60_id_indicator.GMA_lps.irf.(f)...
    + (1-S1_60_id_indicator.aLP_RS.alpha.(f)) .* S1_60_id_indicator.GMA_vars_girf.irf.(f); 
    S1_60_id_indicator.aLP_RS_GMA_girf.mu.(f) = mean(S1_60_id_indicator.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.aLP_RS_GMA_girf.band_lo.(f) = prctile(S1_60_id_indicator.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S1_60_id_indicator.aLP_RS_GMA_girf.band_hi.(f) = prctile(S1_60_id_indicator.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_60_id_indicator.aLP_RS_CVA_girf.irf.(f) = S1_60_id_indicator.aLP_RS.alpha.(f) .* S1_60_id_indicator.CVA_lps.irf.(f)...
    + (1-S1_60_id_indicator.aLP_RS.alpha.(f)) .* S1_60_id_indicator.CVA_vars_girf.irf.(f); 
    S1_60_id_indicator.aLP_RS_CVA_girf.mu.(f) = mean(S1_60_id_indicator.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.aLP_RS_CVA_girf.band_lo.(f) = prctile(S1_60_id_indicator.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S1_60_id_indicator.aLP_RS_CVA_girf.band_hi.(f) = prctile(S1_60_id_indicator.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S1_60_id_indicator.aLP_MSPE_GMA_girf.irf.(f) = S1_60_id_indicator.aLP_MSPE.alpha.(f) .* S1_60_id_indicator.GMA_lps.irf.(f)...
    + (1-S1_60_id_indicator.aLP_MSPE.alpha.(f)) .* S1_60_id_indicator.GMA_vars_girf.irf.(f); 
    S1_60_id_indicator.aLP_MSPE_GMA_girf.mu.(f) = mean(S1_60_id_indicator.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S1_60_id_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S1_60_id_indicator.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S1_60_id_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_60_id_indicator.aLP_MSPE_CVA_girf.irf.(f) = S1_60_id_indicator.aLP_MSPE.alpha.(f) .* S1_60_id_indicator.CVA_lps.irf.(f)...
    + (1-S1_60_id_indicator.aLP_MSPE.alpha.(f)) .* S1_60_id_indicator.CVA_vars_girf.irf.(f); 
    S1_60_id_indicator.aLP_MSPE_CVA_girf.mu.(f) = mean(S1_60_id_indicator.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S1_60_id_indicator.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S1_60_id_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S1_60_id_indicator.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S1_60_id_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S1_60_id_indicator','-append')
%% ------------------------------------------------------- 
% 1-2. the rolling window (60 months) correlation between inflation (yoy) and
% consumption growth (per capita)
% S1_60_pc_indcator
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\indicator\lag12\Bloom_S1_60_pc_indicator.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S1_60_pc_indicator.(names{i}).irf.(f) = irfpaths.(names{i}).S1_60_pc.(f);
        S1_60_pc_indicator.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S1_60_pc.(f),2);
        S1_60_pc_indicator.(names{i}).band_hi.(f) = band_hi.(names{i}).S1_60_pc.(f);
        S1_60_pc_indicator.(names{i}).band_lo.(f) = band_lo.(names{i}).S1_60_pc.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S1_60_pc.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S1_60_pc_indicator.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_60_pc.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S1_60_pc_indicator.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_60_pc.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S1_60_pc_indicator.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S1_60_pc.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S1_60_pc_indicator.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S1_60_pc.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S1_60_pc_indicator.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S1_60_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S1_60_pc_indicator.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_indicator.GMA_lps.irf.(f) = wirf_val;% H x B
    S1_60_pc_indicator.GMA_lps.mu.(f) = mean(S1_60_pc_indicator.GMA_lps.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_indicator.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S1_60_pc_indicator.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_indicator.CVA_lps.irf.(f) = wirf_val;% H x B
    S1_60_pc_indicator.CVA_lps.mu.(f) = mean(S1_60_pc_indicator.CVA_lps.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_indicator.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S1_60_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_60_pc_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_indicator.GMA_vars.irf.(f) = wirf_val;% H x B
    S1_60_pc_indicator.GMA_vars.mu.(f) = mean(S1_60_pc_indicator.GMA_vars.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_indicator.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_60_pc_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_indicator.CVA_vars.irf.(f) = wirf_val;
    S1_60_pc_indicator.CVA_vars.mu.(f) = mean(S1_60_pc_indicator.CVA_vars.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_indicator.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S1_60_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_60_pc_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_indicator.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S1_60_pc_indicator.GMA_vars_girf.mu.(f) = mean(S1_60_pc_indicator.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_indicator.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_60_pc_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_60_pc_indicator.CVA_vars_girf.irf.(f) = wirf_val;
    S1_60_pc_indicator.CVA_vars_girf.mu.(f) = mean(S1_60_pc_indicator.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_60_pc_indicator.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S1_60_pc.state0_vol_empl = tss.S1_60_pc.state0_lps_empl;
    tss.S1_60_pc.state0_ff_empl = tss.S1_60_pc.state0_lps_empl;
    tss.S1_60_pc.state1_vol_empl = tss.S1_60_pc.state1_lps_empl;
    tss.S1_60_pc.state1_ff_empl = tss.S1_60_pc.state1_lps_empl;
    tss.S1_60_pc.state0_vol_ip = tss.S1_60_pc.state0_lps_ip;
    tss.S1_60_pc.state0_ff_ip = tss.S1_60_pc.state0_lps_ip;
    tss.S1_60_pc.state1_vol_ip = tss.S1_60_pc.state1_lps_ip;
    tss.S1_60_pc.state1_ff_ip = tss.S1_60_pc.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S1_60_pc.(f) ./ tss.S1_60_pc.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S1_60_pc.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S1_60_pc.state0_vol_empl = tss.S1_60_pc.state0_vars_empl;
    tss.S1_60_pc.state0_ff_empl = tss.S1_60_pc.state0_vars_empl;
    tss.S1_60_pc.state1_vol_empl = tss.S1_60_pc.state1_vars_empl;
    tss.S1_60_pc.state1_ff_empl = tss.S1_60_pc.state1_vars_empl;
    tss.S1_60_pc.state0_vol_ip = tss.S1_60_pc.state0_vars_ip;
    tss.S1_60_pc.state0_ff_ip = tss.S1_60_pc.state0_vars_ip;
    tss.S1_60_pc.state1_vol_ip = tss.S1_60_pc.state1_vars_ip;
    tss.S1_60_pc.state1_ff_ip = tss.S1_60_pc.state1_vars_ip;
    tss.S1_60_pc.linear_vol_empl = tss.S1_60_pc.linear_vars_empl;
    tss.S1_60_pc.linear_vol_ip = tss.S1_60_pc.linear_vars_ip;
    tss.S1_60_pc.linear_ff_empl = tss.S1_60_pc.linear_vars_empl;
    tss.S1_60_pc.linear_ff_ip = tss.S1_60_pc.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S1_60_pc.(f) ./ tss.S1_60_pc.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S1_60_pc.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S1_60_pc_indicator.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S1_60_pc_indicator.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S1_60_pc_indicator.aLP_RS.alpha.(f) .* S1_60_pc_indicator.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S1_60_pc_indicator.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_60_pc_indicator.GMA_vars.weight.(f)'; % H x 2
    S1_60_pc_indicator.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S1_60_pc_indicator.aLP_RS.alpha.(f) .* S1_60_pc_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_60_pc_indicator.CVA_vars.weight.(f)';
    S1_60_pc_indicator.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S1_60_pc_indicator.aLP_MSPE.alpha.(f) .* S1_60_pc_indicator.GMA_lps.weight.(f)';
    alpha_var = 1 - S1_60_pc_indicator.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_60_pc_indicator.GMA_vars.weight.(f)';
    S1_60_pc_indicator.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S1_60_pc_indicator.aLP_MSPE.alpha.(f) .* S1_60_pc_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_60_pc_indicator.CVA_vars.weight.(f)';
    S1_60_pc_indicator.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S1_60_pc_indicator.aLP_RS_GMA.irf.(f) = S1_60_pc_indicator.aLP_RS.alpha.(f) .* S1_60_pc_indicator.GMA_lps.irf.(f)...
    + (1-S1_60_pc_indicator.aLP_RS.alpha.(f)) .* S1_60_pc_indicator.GMA_vars.irf.(f); 
    S1_60_pc_indicator.aLP_RS_GMA.mu.(f) = mean(S1_60_pc_indicator.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.aLP_RS_GMA.band_lo.(f) = prctile(S1_60_pc_indicator.aLP_RS_GMA.irf.(f),pct_lo,2);
    S1_60_pc_indicator.aLP_RS_GMA.band_hi.(f) = prctile(S1_60_pc_indicator.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S1_60_pc_indicator.aLP_RS_CVA.irf.(f) = S1_60_pc_indicator.aLP_RS.alpha.(f) .* S1_60_pc_indicator.CVA_lps.irf.(f)...
    + (1-S1_60_pc_indicator.aLP_RS.alpha.(f)) .* S1_60_pc_indicator.CVA_vars.irf.(f); 
    S1_60_pc_indicator.aLP_RS_CVA.mu.(f) = mean(S1_60_pc_indicator.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.aLP_RS_CVA.band_lo.(f) = prctile(S1_60_pc_indicator.aLP_RS_CVA.irf.(f),pct_lo,2);
    S1_60_pc_indicator.aLP_RS_CVA.band_hi.(f) = prctile(S1_60_pc_indicator.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S1_60_pc_indicator.aLP_MSPE_GMA.irf.(f) = S1_60_pc_indicator.aLP_MSPE.alpha.(f) .* S1_60_pc_indicator.GMA_lps.irf.(f)...
    + (1-S1_60_pc_indicator.aLP_MSPE.alpha.(f)) .* S1_60_pc_indicator.GMA_vars.irf.(f); 
    S1_60_pc_indicator.aLP_MSPE_GMA.mu.(f) = mean(S1_60_pc_indicator.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.aLP_MSPE_GMA.band_lo.(f) = prctile(S1_60_pc_indicator.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S1_60_pc_indicator.aLP_MSPE_GMA.band_hi.(f) = prctile(S1_60_pc_indicator.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S1_60_pc_indicator.aLP_MSPE_CVA.irf.(f) = S1_60_pc_indicator.aLP_MSPE.alpha.(f) .* S1_60_pc_indicator.CVA_lps.irf.(f)...
    + (1-S1_60_pc_indicator.aLP_MSPE.alpha.(f)) .* S1_60_pc_indicator.CVA_vars.irf.(f); 
    S1_60_pc_indicator.aLP_MSPE_CVA.mu.(f) = mean(S1_60_pc_indicator.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.aLP_MSPE_CVA.band_lo.(f) = prctile(S1_60_pc_indicator.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S1_60_pc_indicator.aLP_MSPE_CVA.band_hi.(f) = prctile(S1_60_pc_indicator.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S1_60_pc_indicator.aLP_RS_GMA_girf.irf.(f) = S1_60_pc_indicator.aLP_RS.alpha.(f) .* S1_60_pc_indicator.GMA_lps.irf.(f)...
    + (1-S1_60_pc_indicator.aLP_RS.alpha.(f)) .* S1_60_pc_indicator.GMA_vars_girf.irf.(f); 
    S1_60_pc_indicator.aLP_RS_GMA_girf.mu.(f) = mean(S1_60_pc_indicator.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.aLP_RS_GMA_girf.band_lo.(f) = prctile(S1_60_pc_indicator.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S1_60_pc_indicator.aLP_RS_GMA_girf.band_hi.(f) = prctile(S1_60_pc_indicator.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_60_pc_indicator.aLP_RS_CVA_girf.irf.(f) = S1_60_pc_indicator.aLP_RS.alpha.(f) .* S1_60_pc_indicator.CVA_lps.irf.(f)...
    + (1-S1_60_pc_indicator.aLP_RS.alpha.(f)) .* S1_60_pc_indicator.CVA_vars_girf.irf.(f); 
    S1_60_pc_indicator.aLP_RS_CVA_girf.mu.(f) = mean(S1_60_pc_indicator.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.aLP_RS_CVA_girf.band_lo.(f) = prctile(S1_60_pc_indicator.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S1_60_pc_indicator.aLP_RS_CVA_girf.band_hi.(f) = prctile(S1_60_pc_indicator.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S1_60_pc_indicator.aLP_MSPE_GMA_girf.irf.(f) = S1_60_pc_indicator.aLP_MSPE.alpha.(f) .* S1_60_pc_indicator.GMA_lps.irf.(f)...
    + (1-S1_60_pc_indicator.aLP_MSPE.alpha.(f)) .* S1_60_pc_indicator.GMA_vars_girf.irf.(f); 
    S1_60_pc_indicator.aLP_MSPE_GMA_girf.mu.(f) = mean(S1_60_pc_indicator.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S1_60_pc_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S1_60_pc_indicator.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S1_60_pc_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_60_pc_indicator.aLP_MSPE_CVA_girf.irf.(f) = S1_60_pc_indicator.aLP_MSPE.alpha.(f) .* S1_60_pc_indicator.CVA_lps.irf.(f)...
    + (1-S1_60_pc_indicator.aLP_MSPE.alpha.(f)) .* S1_60_pc_indicator.CVA_vars_girf.irf.(f); 
    S1_60_pc_indicator.aLP_MSPE_CVA_girf.mu.(f) = mean(S1_60_pc_indicator.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S1_60_pc_indicator.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S1_60_pc_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S1_60_pc_indicator.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S1_60_pc_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S1_60_pc_indicator','-append')
%% ------------------------------------------------------- 
% 1-3. the rolling window (30 months) correlation between inflation (yoy) and
% consumption growth (index)
% S1_30_id_indcator
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\indicator\lag12\Bloom_S1_30_id_indicator.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S1_30_id_indicator.(names{i}).irf.(f) = irfpaths.(names{i}).S1_30_id.(f);
        S1_30_id_indicator.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S1_30_id.(f),2);
        S1_30_id_indicator.(names{i}).band_hi.(f) = band_hi.(names{i}).S1_30_id.(f);
        S1_30_id_indicator.(names{i}).band_lo.(f) = band_lo.(names{i}).S1_30_id.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S1_30_id.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S1_30_id_indicator.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_30_id.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S1_30_id_indicator.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_30_id.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S1_30_id_indicator.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S1_30_id.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S1_30_id_indicator.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S1_30_id.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S1_30_id_indicator.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S1_30_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S1_30_id_indicator.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_indicator.GMA_lps.irf.(f) = wirf_val;% H x B
    S1_30_id_indicator.GMA_lps.mu.(f) = mean(S1_30_id_indicator.GMA_lps.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_indicator.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S1_30_id_indicator.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_indicator.CVA_lps.irf.(f) = wirf_val;% H x B
    S1_30_id_indicator.CVA_lps.mu.(f) = mean(S1_30_id_indicator.CVA_lps.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_indicator.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S1_30_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_30_id_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_indicator.GMA_vars.irf.(f) = wirf_val;% H x B
    S1_30_id_indicator.GMA_vars.mu.(f) = mean(S1_30_id_indicator.GMA_vars.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_indicator.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_30_id_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_indicator.CVA_vars.irf.(f) = wirf_val;
    S1_30_id_indicator.CVA_vars.mu.(f) = mean(S1_30_id_indicator.CVA_vars.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_indicator.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S1_30_id.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_30_id_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_indicator.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S1_30_id_indicator.GMA_vars_girf.mu.(f) = mean(S1_30_id_indicator.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_indicator.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_30_id_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_id_indicator.CVA_vars_girf.irf.(f) = wirf_val;
    S1_30_id_indicator.CVA_vars_girf.mu.(f) = mean(S1_30_id_indicator.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_id_indicator.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S1_30_id.state0_vol_empl = tss.S1_30_id.state0_lps_empl;
    tss.S1_30_id.state0_ff_empl = tss.S1_30_id.state0_lps_empl;
    tss.S1_30_id.state1_vol_empl = tss.S1_30_id.state1_lps_empl;
    tss.S1_30_id.state1_ff_empl = tss.S1_30_id.state1_lps_empl;
    tss.S1_30_id.state0_vol_ip = tss.S1_30_id.state0_lps_ip;
    tss.S1_30_id.state0_ff_ip = tss.S1_30_id.state0_lps_ip;
    tss.S1_30_id.state1_vol_ip = tss.S1_30_id.state1_lps_ip;
    tss.S1_30_id.state1_ff_ip = tss.S1_30_id.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S1_30_id.(f) ./ tss.S1_30_id.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S1_30_id.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S1_30_id.state0_vol_empl = tss.S1_30_id.state0_vars_empl;
    tss.S1_30_id.state0_ff_empl = tss.S1_30_id.state0_vars_empl;
    tss.S1_30_id.state1_vol_empl = tss.S1_30_id.state1_vars_empl;
    tss.S1_30_id.state1_ff_empl = tss.S1_30_id.state1_vars_empl;
    tss.S1_30_id.state0_vol_ip = tss.S1_30_id.state0_vars_ip;
    tss.S1_30_id.state0_ff_ip = tss.S1_30_id.state0_vars_ip;
    tss.S1_30_id.state1_vol_ip = tss.S1_30_id.state1_vars_ip;
    tss.S1_30_id.state1_ff_ip = tss.S1_30_id.state1_vars_ip;
    tss.S1_30_id.linear_vol_empl = tss.S1_30_id.linear_vars_empl;
    tss.S1_30_id.linear_vol_ip = tss.S1_30_id.linear_vars_ip;
    tss.S1_30_id.linear_ff_empl = tss.S1_30_id.linear_vars_empl;
    tss.S1_30_id.linear_ff_ip = tss.S1_30_id.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S1_30_id.(f) ./ tss.S1_30_id.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S1_30_id.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S1_30_id_indicator.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S1_30_id_indicator.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S1_30_id_indicator.aLP_RS.alpha.(f) .* S1_30_id_indicator.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S1_30_id_indicator.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_30_id_indicator.GMA_vars.weight.(f)'; % H x 2
    S1_30_id_indicator.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S1_30_id_indicator.aLP_RS.alpha.(f) .* S1_30_id_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_30_id_indicator.CVA_vars.weight.(f)';
    S1_30_id_indicator.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S1_30_id_indicator.aLP_MSPE.alpha.(f) .* S1_30_id_indicator.GMA_lps.weight.(f)';
    alpha_var = 1 - S1_30_id_indicator.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_30_id_indicator.GMA_vars.weight.(f)';
    S1_30_id_indicator.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S1_30_id_indicator.aLP_MSPE.alpha.(f) .* S1_30_id_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_30_id_indicator.CVA_vars.weight.(f)';
    S1_30_id_indicator.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S1_30_id_indicator.aLP_RS_GMA.irf.(f) = S1_30_id_indicator.aLP_RS.alpha.(f) .* S1_30_id_indicator.GMA_lps.irf.(f)...
    + (1-S1_30_id_indicator.aLP_RS.alpha.(f)) .* S1_30_id_indicator.GMA_vars.irf.(f); 
    S1_30_id_indicator.aLP_RS_GMA.mu.(f) = mean(S1_30_id_indicator.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.aLP_RS_GMA.band_lo.(f) = prctile(S1_30_id_indicator.aLP_RS_GMA.irf.(f),pct_lo,2);
    S1_30_id_indicator.aLP_RS_GMA.band_hi.(f) = prctile(S1_30_id_indicator.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S1_30_id_indicator.aLP_RS_CVA.irf.(f) = S1_30_id_indicator.aLP_RS.alpha.(f) .* S1_30_id_indicator.CVA_lps.irf.(f)...
    + (1-S1_30_id_indicator.aLP_RS.alpha.(f)) .* S1_30_id_indicator.CVA_vars.irf.(f); 
    S1_30_id_indicator.aLP_RS_CVA.mu.(f) = mean(S1_30_id_indicator.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.aLP_RS_CVA.band_lo.(f) = prctile(S1_30_id_indicator.aLP_RS_CVA.irf.(f),pct_lo,2);
    S1_30_id_indicator.aLP_RS_CVA.band_hi.(f) = prctile(S1_30_id_indicator.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S1_30_id_indicator.aLP_MSPE_GMA.irf.(f) = S1_30_id_indicator.aLP_MSPE.alpha.(f) .* S1_30_id_indicator.GMA_lps.irf.(f)...
    + (1-S1_30_id_indicator.aLP_MSPE.alpha.(f)) .* S1_30_id_indicator.GMA_vars.irf.(f); 
    S1_30_id_indicator.aLP_MSPE_GMA.mu.(f) = mean(S1_30_id_indicator.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.aLP_MSPE_GMA.band_lo.(f) = prctile(S1_30_id_indicator.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S1_30_id_indicator.aLP_MSPE_GMA.band_hi.(f) = prctile(S1_30_id_indicator.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S1_30_id_indicator.aLP_MSPE_CVA.irf.(f) = S1_30_id_indicator.aLP_MSPE.alpha.(f) .* S1_30_id_indicator.CVA_lps.irf.(f)...
    + (1-S1_30_id_indicator.aLP_MSPE.alpha.(f)) .* S1_30_id_indicator.CVA_vars.irf.(f); 
    S1_30_id_indicator.aLP_MSPE_CVA.mu.(f) = mean(S1_30_id_indicator.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.aLP_MSPE_CVA.band_lo.(f) = prctile(S1_30_id_indicator.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S1_30_id_indicator.aLP_MSPE_CVA.band_hi.(f) = prctile(S1_30_id_indicator.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S1_30_id_indicator.aLP_RS_GMA_girf.irf.(f) = S1_30_id_indicator.aLP_RS.alpha.(f) .* S1_30_id_indicator.GMA_lps.irf.(f)...
    + (1-S1_30_id_indicator.aLP_RS.alpha.(f)) .* S1_30_id_indicator.GMA_vars_girf.irf.(f); 
    S1_30_id_indicator.aLP_RS_GMA_girf.mu.(f) = mean(S1_30_id_indicator.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.aLP_RS_GMA_girf.band_lo.(f) = prctile(S1_30_id_indicator.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S1_30_id_indicator.aLP_RS_GMA_girf.band_hi.(f) = prctile(S1_30_id_indicator.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_30_id_indicator.aLP_RS_CVA_girf.irf.(f) = S1_30_id_indicator.aLP_RS.alpha.(f) .* S1_30_id_indicator.CVA_lps.irf.(f)...
    + (1-S1_30_id_indicator.aLP_RS.alpha.(f)) .* S1_30_id_indicator.CVA_vars_girf.irf.(f); 
    S1_30_id_indicator.aLP_RS_CVA_girf.mu.(f) = mean(S1_30_id_indicator.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.aLP_RS_CVA_girf.band_lo.(f) = prctile(S1_30_id_indicator.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S1_30_id_indicator.aLP_RS_CVA_girf.band_hi.(f) = prctile(S1_30_id_indicator.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S1_30_id_indicator.aLP_MSPE_GMA_girf.irf.(f) = S1_30_id_indicator.aLP_MSPE.alpha.(f) .* S1_30_id_indicator.GMA_lps.irf.(f)...
    + (1-S1_30_id_indicator.aLP_MSPE.alpha.(f)) .* S1_30_id_indicator.GMA_vars_girf.irf.(f); 
    S1_30_id_indicator.aLP_MSPE_GMA_girf.mu.(f) = mean(S1_30_id_indicator.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S1_30_id_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S1_30_id_indicator.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S1_30_id_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_30_id_indicator.aLP_MSPE_CVA_girf.irf.(f) = S1_30_id_indicator.aLP_MSPE.alpha.(f) .* S1_30_id_indicator.CVA_lps.irf.(f)...
    + (1-S1_30_id_indicator.aLP_MSPE.alpha.(f)) .* S1_30_id_indicator.CVA_vars_girf.irf.(f); 
    S1_30_id_indicator.aLP_MSPE_CVA_girf.mu.(f) = mean(S1_30_id_indicator.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S1_30_id_indicator.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S1_30_id_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S1_30_id_indicator.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S1_30_id_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S1_30_id_indicator','-append')
%% ------------------------------------------------------- 
% 1-4. the rolling window (30 months) correlation between inflation (yoy) and
% consumption growth (per capita)
% S1_30_pc_indcator
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\indicator\lag12\Bloom_S1_30_pc_indicator.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S1_30_pc_indicator.(names{i}).irf.(f) = irfpaths.(names{i}).S1_30_pc.(f);
        S1_30_pc_indicator.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S1_30_pc.(f),2);
        S1_30_pc_indicator.(names{i}).band_hi.(f) = band_hi.(names{i}).S1_30_pc.(f);
        S1_30_pc_indicator.(names{i}).band_lo.(f) = band_lo.(names{i}).S1_30_pc.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S1_30_pc.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S1_30_pc_indicator.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_30_pc.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S1_30_pc_indicator.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S1_30_pc.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S1_30_pc_indicator.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S1_30_pc.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S1_30_pc_indicator.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S1_30_pc.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S1_30_pc_indicator.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S1_30_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S1_30_pc_indicator.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_indicator.GMA_lps.irf.(f) = wirf_val;% H x B
    S1_30_pc_indicator.GMA_lps.mu.(f) = mean(S1_30_pc_indicator.GMA_lps.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_indicator.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S1_30_pc_indicator.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_indicator.CVA_lps.irf.(f) = wirf_val;% H x B
    S1_30_pc_indicator.CVA_lps.mu.(f) = mean(S1_30_pc_indicator.CVA_lps.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_indicator.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S1_30_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_30_pc_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_indicator.GMA_vars.irf.(f) = wirf_val;% H x B
    S1_30_pc_indicator.GMA_vars.mu.(f) = mean(S1_30_pc_indicator.GMA_vars.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_indicator.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_30_pc_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_indicator.CVA_vars.irf.(f) = wirf_val;
    S1_30_pc_indicator.CVA_vars.mu.(f) = mean(S1_30_pc_indicator.CVA_vars.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_indicator.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S1_30_pc.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S1_30_pc_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_indicator.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S1_30_pc_indicator.GMA_vars_girf.mu.(f) = mean(S1_30_pc_indicator.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_indicator.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S1_30_pc_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S1_30_pc_indicator.CVA_vars_girf.irf.(f) = wirf_val;
    S1_30_pc_indicator.CVA_vars_girf.mu.(f) = mean(S1_30_pc_indicator.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S1_30_pc_indicator.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S1_30_pc.state0_vol_empl = tss.S1_30_pc.state0_lps_empl;
    tss.S1_30_pc.state0_ff_empl = tss.S1_30_pc.state0_lps_empl;
    tss.S1_30_pc.state1_vol_empl = tss.S1_30_pc.state1_lps_empl;
    tss.S1_30_pc.state1_ff_empl = tss.S1_30_pc.state1_lps_empl;
    tss.S1_30_pc.state0_vol_ip = tss.S1_30_pc.state0_lps_ip;
    tss.S1_30_pc.state0_ff_ip = tss.S1_30_pc.state0_lps_ip;
    tss.S1_30_pc.state1_vol_ip = tss.S1_30_pc.state1_lps_ip;
    tss.S1_30_pc.state1_ff_ip = tss.S1_30_pc.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S1_30_pc.(f) ./ tss.S1_30_pc.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S1_30_pc.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S1_30_pc.state0_vol_empl = tss.S1_30_pc.state0_vars_empl;
    tss.S1_30_pc.state0_ff_empl = tss.S1_30_pc.state0_vars_empl;
    tss.S1_30_pc.state1_vol_empl = tss.S1_30_pc.state1_vars_empl;
    tss.S1_30_pc.state1_ff_empl = tss.S1_30_pc.state1_vars_empl;
    tss.S1_30_pc.state0_vol_ip = tss.S1_30_pc.state0_vars_ip;
    tss.S1_30_pc.state0_ff_ip = tss.S1_30_pc.state0_vars_ip;
    tss.S1_30_pc.state1_vol_ip = tss.S1_30_pc.state1_vars_ip;
    tss.S1_30_pc.state1_ff_ip = tss.S1_30_pc.state1_vars_ip;
    tss.S1_30_pc.linear_vol_empl = tss.S1_30_pc.linear_vars_empl;
    tss.S1_30_pc.linear_vol_ip = tss.S1_30_pc.linear_vars_ip;
    tss.S1_30_pc.linear_ff_empl = tss.S1_30_pc.linear_vars_empl;
    tss.S1_30_pc.linear_ff_ip = tss.S1_30_pc.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S1_30_pc.(f) ./ tss.S1_30_pc.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S1_30_pc.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S1_30_pc_indicator.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S1_30_pc_indicator.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S1_30_pc_indicator.aLP_RS.alpha.(f) .* S1_30_pc_indicator.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S1_30_pc_indicator.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_30_pc_indicator.GMA_vars.weight.(f)'; % H x 2
    S1_30_pc_indicator.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S1_30_pc_indicator.aLP_RS.alpha.(f) .* S1_30_pc_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_30_pc_indicator.CVA_vars.weight.(f)';
    S1_30_pc_indicator.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S1_30_pc_indicator.aLP_MSPE.alpha.(f) .* S1_30_pc_indicator.GMA_lps.weight.(f)';
    alpha_var = 1 - S1_30_pc_indicator.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S1_30_pc_indicator.GMA_vars.weight.(f)';
    S1_30_pc_indicator.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S1_30_pc_indicator.aLP_MSPE.alpha.(f) .* S1_30_pc_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S1_30_pc_indicator.CVA_vars.weight.(f)';
    S1_30_pc_indicator.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S1_30_pc_indicator.aLP_RS_GMA.irf.(f) = S1_30_pc_indicator.aLP_RS.alpha.(f) .* S1_30_pc_indicator.GMA_lps.irf.(f)...
    + (1-S1_30_pc_indicator.aLP_RS.alpha.(f)) .* S1_30_pc_indicator.GMA_vars.irf.(f); 
    S1_30_pc_indicator.aLP_RS_GMA.mu.(f) = mean(S1_30_pc_indicator.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.aLP_RS_GMA.band_lo.(f) = prctile(S1_30_pc_indicator.aLP_RS_GMA.irf.(f),pct_lo,2);
    S1_30_pc_indicator.aLP_RS_GMA.band_hi.(f) = prctile(S1_30_pc_indicator.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S1_30_pc_indicator.aLP_RS_CVA.irf.(f) = S1_30_pc_indicator.aLP_RS.alpha.(f) .* S1_30_pc_indicator.CVA_lps.irf.(f)...
    + (1-S1_30_pc_indicator.aLP_RS.alpha.(f)) .* S1_30_pc_indicator.CVA_vars.irf.(f); 
    S1_30_pc_indicator.aLP_RS_CVA.mu.(f) = mean(S1_30_pc_indicator.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.aLP_RS_CVA.band_lo.(f) = prctile(S1_30_pc_indicator.aLP_RS_CVA.irf.(f),pct_lo,2);
    S1_30_pc_indicator.aLP_RS_CVA.band_hi.(f) = prctile(S1_30_pc_indicator.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S1_30_pc_indicator.aLP_MSPE_GMA.irf.(f) = S1_30_pc_indicator.aLP_MSPE.alpha.(f) .* S1_30_pc_indicator.GMA_lps.irf.(f)...
    + (1-S1_30_pc_indicator.aLP_MSPE.alpha.(f)) .* S1_30_pc_indicator.GMA_vars.irf.(f); 
    S1_30_pc_indicator.aLP_MSPE_GMA.mu.(f) = mean(S1_30_pc_indicator.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.aLP_MSPE_GMA.band_lo.(f) = prctile(S1_30_pc_indicator.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S1_30_pc_indicator.aLP_MSPE_GMA.band_hi.(f) = prctile(S1_30_pc_indicator.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S1_30_pc_indicator.aLP_MSPE_CVA.irf.(f) = S1_30_pc_indicator.aLP_MSPE.alpha.(f) .* S1_30_pc_indicator.CVA_lps.irf.(f)...
    + (1-S1_30_pc_indicator.aLP_MSPE.alpha.(f)) .* S1_30_pc_indicator.CVA_vars.irf.(f); 
    S1_30_pc_indicator.aLP_MSPE_CVA.mu.(f) = mean(S1_30_pc_indicator.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.aLP_MSPE_CVA.band_lo.(f) = prctile(S1_30_pc_indicator.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S1_30_pc_indicator.aLP_MSPE_CVA.band_hi.(f) = prctile(S1_30_pc_indicator.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S1_30_pc_indicator.aLP_RS_GMA_girf.irf.(f) = S1_30_pc_indicator.aLP_RS.alpha.(f) .* S1_30_pc_indicator.GMA_lps.irf.(f)...
    + (1-S1_30_pc_indicator.aLP_RS.alpha.(f)) .* S1_30_pc_indicator.GMA_vars_girf.irf.(f); 
    S1_30_pc_indicator.aLP_RS_GMA_girf.mu.(f) = mean(S1_30_pc_indicator.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.aLP_RS_GMA_girf.band_lo.(f) = prctile(S1_30_pc_indicator.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S1_30_pc_indicator.aLP_RS_GMA_girf.band_hi.(f) = prctile(S1_30_pc_indicator.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_30_pc_indicator.aLP_RS_CVA_girf.irf.(f) = S1_30_pc_indicator.aLP_RS.alpha.(f) .* S1_30_pc_indicator.CVA_lps.irf.(f)...
    + (1-S1_30_pc_indicator.aLP_RS.alpha.(f)) .* S1_30_pc_indicator.CVA_vars_girf.irf.(f); 
    S1_30_pc_indicator.aLP_RS_CVA_girf.mu.(f) = mean(S1_30_pc_indicator.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.aLP_RS_CVA_girf.band_lo.(f) = prctile(S1_30_pc_indicator.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S1_30_pc_indicator.aLP_RS_CVA_girf.band_hi.(f) = prctile(S1_30_pc_indicator.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S1_30_pc_indicator.aLP_MSPE_GMA_girf.irf.(f) = S1_30_pc_indicator.aLP_MSPE.alpha.(f) .* S1_30_pc_indicator.GMA_lps.irf.(f)...
    + (1-S1_30_pc_indicator.aLP_MSPE.alpha.(f)) .* S1_30_pc_indicator.GMA_vars_girf.irf.(f); 
    S1_30_pc_indicator.aLP_MSPE_GMA_girf.mu.(f) = mean(S1_30_pc_indicator.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S1_30_pc_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S1_30_pc_indicator.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S1_30_pc_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S1_30_pc_indicator.aLP_MSPE_CVA_girf.irf.(f) = S1_30_pc_indicator.aLP_MSPE.alpha.(f) .* S1_30_pc_indicator.CVA_lps.irf.(f)...
    + (1-S1_30_pc_indicator.aLP_MSPE.alpha.(f)) .* S1_30_pc_indicator.CVA_vars_girf.irf.(f); 
    S1_30_pc_indicator.aLP_MSPE_CVA_girf.mu.(f) = mean(S1_30_pc_indicator.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S1_30_pc_indicator.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S1_30_pc_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S1_30_pc_indicator.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S1_30_pc_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S1_30_pc_indicator','-append')
%% -------------------------------------------------------
% 2. Anchored vs. Unanchored expectations
% 2-1. Michigan Household Survey IQR from jan 1987
% S2_78_hh_indcator
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\indicator\lag12\Bloom_S2_78_hh_indicator.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S2_78_hh_indicator.(names{i}).irf.(f) = irfpaths.(names{i}).S2_78_hh.(f);
        S2_78_hh_indicator.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S2_78_hh.(f),2);
        S2_78_hh_indicator.(names{i}).band_hi.(f) = band_hi.(names{i}).S2_78_hh.(f);
        S2_78_hh_indicator.(names{i}).band_lo.(f) = band_lo.(names{i}).S2_78_hh.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S2_78_hh.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S2_78_hh_indicator.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_78_hh.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S2_78_hh_indicator.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_78_hh.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S2_78_hh_indicator.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S2_78_hh.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S2_78_hh_indicator.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S2_78_hh.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S2_78_hh_indicator.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S2_78_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S2_78_hh_indicator.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_indicator.GMA_lps.irf.(f) = wirf_val;% H x B
    S2_78_hh_indicator.GMA_lps.mu.(f) = mean(S2_78_hh_indicator.GMA_lps.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_indicator.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S2_78_hh_indicator.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_indicator.CVA_lps.irf.(f) = wirf_val;% H x B
    S2_78_hh_indicator.CVA_lps.mu.(f) = mean(S2_78_hh_indicator.CVA_lps.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_indicator.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S2_78_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_78_hh_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_indicator.GMA_vars.irf.(f) = wirf_val;% H x B
    S2_78_hh_indicator.GMA_vars.mu.(f) = mean(S2_78_hh_indicator.GMA_vars.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_indicator.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_78_hh_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_indicator.CVA_vars.irf.(f) = wirf_val;
    S2_78_hh_indicator.CVA_vars.mu.(f) = mean(S2_78_hh_indicator.CVA_vars.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_indicator.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S2_78_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_78_hh_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_indicator.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S2_78_hh_indicator.GMA_vars_girf.mu.(f) = mean(S2_78_hh_indicator.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_indicator.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_78_hh_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_78_hh_indicator.CVA_vars_girf.irf.(f) = wirf_val;
    S2_78_hh_indicator.CVA_vars_girf.mu.(f) = mean(S2_78_hh_indicator.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_78_hh_indicator.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S2_78_hh.state0_vol_empl = tss.S2_78_hh.state0_lps_empl;
    tss.S2_78_hh.state0_ff_empl = tss.S2_78_hh.state0_lps_empl;
    tss.S2_78_hh.state1_vol_empl = tss.S2_78_hh.state1_lps_empl;
    tss.S2_78_hh.state1_ff_empl = tss.S2_78_hh.state1_lps_empl;
    tss.S2_78_hh.state0_vol_ip = tss.S2_78_hh.state0_lps_ip;
    tss.S2_78_hh.state0_ff_ip = tss.S2_78_hh.state0_lps_ip;
    tss.S2_78_hh.state1_vol_ip = tss.S2_78_hh.state1_lps_ip;
    tss.S2_78_hh.state1_ff_ip = tss.S2_78_hh.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S2_78_hh.(f) ./ tss.S2_78_hh.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S2_78_hh.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S2_78_hh.state0_vol_empl = tss.S2_78_hh.state0_vars_empl;
    tss.S2_78_hh.state0_ff_empl = tss.S2_78_hh.state0_vars_empl;
    tss.S2_78_hh.state1_vol_empl = tss.S2_78_hh.state1_vars_empl;
    tss.S2_78_hh.state1_ff_empl = tss.S2_78_hh.state1_vars_empl;
    tss.S2_78_hh.state0_vol_ip = tss.S2_78_hh.state0_vars_ip;
    tss.S2_78_hh.state0_ff_ip = tss.S2_78_hh.state0_vars_ip;
    tss.S2_78_hh.state1_vol_ip = tss.S2_78_hh.state1_vars_ip;
    tss.S2_78_hh.state1_ff_ip = tss.S2_78_hh.state1_vars_ip;
    tss.S2_78_hh.linear_vol_empl = tss.S2_78_hh.linear_vars_empl;
    tss.S2_78_hh.linear_vol_ip = tss.S2_78_hh.linear_vars_ip;
    tss.S2_78_hh.linear_ff_empl = tss.S2_78_hh.linear_vars_empl;
    tss.S2_78_hh.linear_ff_ip = tss.S2_78_hh.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S2_78_hh.(f) ./ tss.S2_78_hh.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S2_78_hh.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S2_78_hh_indicator.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S2_78_hh_indicator.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S2_78_hh_indicator.aLP_RS.alpha.(f) .* S2_78_hh_indicator.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S2_78_hh_indicator.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_78_hh_indicator.GMA_vars.weight.(f)'; % H x 2
    S2_78_hh_indicator.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S2_78_hh_indicator.aLP_RS.alpha.(f) .* S2_78_hh_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_78_hh_indicator.CVA_vars.weight.(f)';
    S2_78_hh_indicator.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S2_78_hh_indicator.aLP_MSPE.alpha.(f) .* S2_78_hh_indicator.GMA_lps.weight.(f)';
    alpha_var = 1 - S2_78_hh_indicator.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_78_hh_indicator.GMA_vars.weight.(f)';
    S2_78_hh_indicator.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S2_78_hh_indicator.aLP_MSPE.alpha.(f) .* S2_78_hh_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_78_hh_indicator.CVA_vars.weight.(f)';
    S2_78_hh_indicator.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S2_78_hh_indicator.aLP_RS_GMA.irf.(f) = S2_78_hh_indicator.aLP_RS.alpha.(f) .* S2_78_hh_indicator.GMA_lps.irf.(f)...
    + (1-S2_78_hh_indicator.aLP_RS.alpha.(f)) .* S2_78_hh_indicator.GMA_vars.irf.(f); 
    S2_78_hh_indicator.aLP_RS_GMA.mu.(f) = mean(S2_78_hh_indicator.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.aLP_RS_GMA.band_lo.(f) = prctile(S2_78_hh_indicator.aLP_RS_GMA.irf.(f),pct_lo,2);
    S2_78_hh_indicator.aLP_RS_GMA.band_hi.(f) = prctile(S2_78_hh_indicator.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S2_78_hh_indicator.aLP_RS_CVA.irf.(f) = S2_78_hh_indicator.aLP_RS.alpha.(f) .* S2_78_hh_indicator.CVA_lps.irf.(f)...
    + (1-S2_78_hh_indicator.aLP_RS.alpha.(f)) .* S2_78_hh_indicator.CVA_vars.irf.(f); 
    S2_78_hh_indicator.aLP_RS_CVA.mu.(f) = mean(S2_78_hh_indicator.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.aLP_RS_CVA.band_lo.(f) = prctile(S2_78_hh_indicator.aLP_RS_CVA.irf.(f),pct_lo,2);
    S2_78_hh_indicator.aLP_RS_CVA.band_hi.(f) = prctile(S2_78_hh_indicator.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S2_78_hh_indicator.aLP_MSPE_GMA.irf.(f) = S2_78_hh_indicator.aLP_MSPE.alpha.(f) .* S2_78_hh_indicator.GMA_lps.irf.(f)...
    + (1-S2_78_hh_indicator.aLP_MSPE.alpha.(f)) .* S2_78_hh_indicator.GMA_vars.irf.(f); 
    S2_78_hh_indicator.aLP_MSPE_GMA.mu.(f) = mean(S2_78_hh_indicator.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.aLP_MSPE_GMA.band_lo.(f) = prctile(S2_78_hh_indicator.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S2_78_hh_indicator.aLP_MSPE_GMA.band_hi.(f) = prctile(S2_78_hh_indicator.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S2_78_hh_indicator.aLP_MSPE_CVA.irf.(f) = S2_78_hh_indicator.aLP_MSPE.alpha.(f) .* S2_78_hh_indicator.CVA_lps.irf.(f)...
    + (1-S2_78_hh_indicator.aLP_MSPE.alpha.(f)) .* S2_78_hh_indicator.CVA_vars.irf.(f); 
    S2_78_hh_indicator.aLP_MSPE_CVA.mu.(f) = mean(S2_78_hh_indicator.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.aLP_MSPE_CVA.band_lo.(f) = prctile(S2_78_hh_indicator.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S2_78_hh_indicator.aLP_MSPE_CVA.band_hi.(f) = prctile(S2_78_hh_indicator.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S2_78_hh_indicator.aLP_RS_GMA_girf.irf.(f) = S2_78_hh_indicator.aLP_RS.alpha.(f) .* S2_78_hh_indicator.GMA_lps.irf.(f)...
    + (1-S2_78_hh_indicator.aLP_RS.alpha.(f)) .* S2_78_hh_indicator.GMA_vars_girf.irf.(f); 
    S2_78_hh_indicator.aLP_RS_GMA_girf.mu.(f) = mean(S2_78_hh_indicator.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.aLP_RS_GMA_girf.band_lo.(f) = prctile(S2_78_hh_indicator.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S2_78_hh_indicator.aLP_RS_GMA_girf.band_hi.(f) = prctile(S2_78_hh_indicator.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_78_hh_indicator.aLP_RS_CVA_girf.irf.(f) = S2_78_hh_indicator.aLP_RS.alpha.(f) .* S2_78_hh_indicator.CVA_lps.irf.(f)...
    + (1-S2_78_hh_indicator.aLP_RS.alpha.(f)) .* S2_78_hh_indicator.CVA_vars_girf.irf.(f); 
    S2_78_hh_indicator.aLP_RS_CVA_girf.mu.(f) = mean(S2_78_hh_indicator.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.aLP_RS_CVA_girf.band_lo.(f) = prctile(S2_78_hh_indicator.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S2_78_hh_indicator.aLP_RS_CVA_girf.band_hi.(f) = prctile(S2_78_hh_indicator.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S2_78_hh_indicator.aLP_MSPE_GMA_girf.irf.(f) = S2_78_hh_indicator.aLP_MSPE.alpha.(f) .* S2_78_hh_indicator.GMA_lps.irf.(f)...
    + (1-S2_78_hh_indicator.aLP_MSPE.alpha.(f)) .* S2_78_hh_indicator.GMA_vars_girf.irf.(f); 
    S2_78_hh_indicator.aLP_MSPE_GMA_girf.mu.(f) = mean(S2_78_hh_indicator.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S2_78_hh_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S2_78_hh_indicator.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S2_78_hh_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_78_hh_indicator.aLP_MSPE_CVA_girf.irf.(f) = S2_78_hh_indicator.aLP_MSPE.alpha.(f) .* S2_78_hh_indicator.CVA_lps.irf.(f)...
    + (1-S2_78_hh_indicator.aLP_MSPE.alpha.(f)) .* S2_78_hh_indicator.CVA_vars_girf.irf.(f); 
    S2_78_hh_indicator.aLP_MSPE_CVA_girf.mu.(f) = mean(S2_78_hh_indicator.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S2_78_hh_indicator.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S2_78_hh_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S2_78_hh_indicator.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S2_78_hh_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end
%%
save("Combine_Results_Emp.mat",'S2_78_hh_indicator','-append')
%% ------------------------------------------------------
% 2-2. Michigan Household Survey IQR from 1981
% S2_81_hh_indcator: no CVA_lps
% -------------------------------------------------------
    
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\indicator\lag12\Bloom_S2_81_hh_indicator.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S2_81_hh_indicator.(names{i}).irf.(f) = irfpaths.(names{i}).S2_81_hh.(f);
        S2_81_hh_indicator.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S2_81_hh.(f),2);
        S2_81_hh_indicator.(names{i}).band_hi.(f) = band_hi.(names{i}).S2_81_hh.(f);
        S2_81_hh_indicator.(names{i}).band_lo.(f) = band_lo.(names{i}).S2_81_hh.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S2_81_hh.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S2_81_hh_indicator.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_81_hh.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S2_81_hh_indicator.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_81_hh.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S2_81_hh_indicator.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars

% CVA -----------------------------------------------
% settings: S2_81_hh_indicator has no CVA weight for regimes 
for j = 1:length(fields)
    f = fields{j};
    CVA_w_lps.S2_81_hh.(f) = NaN(settings.est.IRF_hor, length(lp_names));
end
%
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S2_81_hh.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S2_81_hh_indicator.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S2_81_hh.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S2_81_hh_indicator.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S2_81_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S2_81_hh_indicator.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_indicator.GMA_lps.irf.(f) = wirf_val;% H x B
    S2_81_hh_indicator.GMA_lps.mu.(f) = mean(S2_81_hh_indicator.GMA_lps.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_indicator.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S2_81_hh_indicator.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_indicator.CVA_lps.irf.(f) = wirf_val;% H x B
    S2_81_hh_indicator.CVA_lps.mu.(f) = mean(S2_81_hh_indicator.CVA_lps.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_indicator.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S2_81_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_81_hh_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_indicator.GMA_vars.irf.(f) = wirf_val;% H x B
    S2_81_hh_indicator.GMA_vars.mu.(f) = mean(S2_81_hh_indicator.GMA_vars.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_indicator.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_81_hh_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_indicator.CVA_vars.irf.(f) = wirf_val;
    S2_81_hh_indicator.CVA_vars.mu.(f) = mean(S2_81_hh_indicator.CVA_vars.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_indicator.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S2_81_hh.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_81_hh_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_indicator.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S2_81_hh_indicator.GMA_vars_girf.mu.(f) = mean(S2_81_hh_indicator.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_indicator.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_81_hh_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_hh_indicator.CVA_vars_girf.irf.(f) = wirf_val;
    S2_81_hh_indicator.CVA_vars_girf.mu.(f) = mean(S2_81_hh_indicator.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_hh_indicator.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S2_81_hh.state0_vol_empl = tss.S2_81_hh.state0_lps_empl;
    tss.S2_81_hh.state0_ff_empl = tss.S2_81_hh.state0_lps_empl;
    tss.S2_81_hh.state1_vol_empl = tss.S2_81_hh.state1_lps_empl;
    tss.S2_81_hh.state1_ff_empl = tss.S2_81_hh.state1_lps_empl;
    tss.S2_81_hh.state0_vol_ip = tss.S2_81_hh.state0_lps_ip;
    tss.S2_81_hh.state0_ff_ip = tss.S2_81_hh.state0_lps_ip;
    tss.S2_81_hh.state1_vol_ip = tss.S2_81_hh.state1_lps_ip;
    tss.S2_81_hh.state1_ff_ip = tss.S2_81_hh.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S2_81_hh.(f) ./ tss.S2_81_hh.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S2_81_hh.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S2_81_hh.state0_vol_empl = tss.S2_81_hh.state0_vars_empl;
    tss.S2_81_hh.state0_ff_empl = tss.S2_81_hh.state0_vars_empl;
    tss.S2_81_hh.state1_vol_empl = tss.S2_81_hh.state1_vars_empl;
    tss.S2_81_hh.state1_ff_empl = tss.S2_81_hh.state1_vars_empl;
    tss.S2_81_hh.state0_vol_ip = tss.S2_81_hh.state0_vars_ip;
    tss.S2_81_hh.state0_ff_ip = tss.S2_81_hh.state0_vars_ip;
    tss.S2_81_hh.state1_vol_ip = tss.S2_81_hh.state1_vars_ip;
    tss.S2_81_hh.state1_ff_ip = tss.S2_81_hh.state1_vars_ip;
    tss.S2_81_hh.linear_vol_empl = tss.S2_81_hh.linear_vars_empl;
    tss.S2_81_hh.linear_vol_ip = tss.S2_81_hh.linear_vars_ip;
    tss.S2_81_hh.linear_ff_empl = tss.S2_81_hh.linear_vars_empl;
    tss.S2_81_hh.linear_ff_ip = tss.S2_81_hh.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S2_81_hh.(f) ./ tss.S2_81_hh.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S2_81_hh.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S2_81_hh_indicator.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S2_81_hh_indicator.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S2_81_hh_indicator.aLP_RS.alpha.(f) .* S2_81_hh_indicator.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S2_81_hh_indicator.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_81_hh_indicator.GMA_vars.weight.(f)'; % H x 2
    S2_81_hh_indicator.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S2_81_hh_indicator.aLP_RS.alpha.(f) .* S2_81_hh_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_81_hh_indicator.CVA_vars.weight.(f)';
    S2_81_hh_indicator.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S2_81_hh_indicator.aLP_MSPE.alpha.(f) .* S2_81_hh_indicator.GMA_lps.weight.(f)';
    alpha_var = 1 - S2_81_hh_indicator.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_81_hh_indicator.GMA_vars.weight.(f)';
    S2_81_hh_indicator.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S2_81_hh_indicator.aLP_MSPE.alpha.(f) .* S2_81_hh_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_81_hh_indicator.CVA_vars.weight.(f)';
    S2_81_hh_indicator.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S2_81_hh_indicator.aLP_RS_GMA.irf.(f) = S2_81_hh_indicator.aLP_RS.alpha.(f) .* S2_81_hh_indicator.GMA_lps.irf.(f)...
    + (1-S2_81_hh_indicator.aLP_RS.alpha.(f)) .* S2_81_hh_indicator.GMA_vars.irf.(f); 
    S2_81_hh_indicator.aLP_RS_GMA.mu.(f) = mean(S2_81_hh_indicator.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.aLP_RS_GMA.band_lo.(f) = prctile(S2_81_hh_indicator.aLP_RS_GMA.irf.(f),pct_lo,2);
    S2_81_hh_indicator.aLP_RS_GMA.band_hi.(f) = prctile(S2_81_hh_indicator.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S2_81_hh_indicator.aLP_RS_CVA.irf.(f) = S2_81_hh_indicator.aLP_RS.alpha.(f) .* S2_81_hh_indicator.CVA_lps.irf.(f)...
    + (1-S2_81_hh_indicator.aLP_RS.alpha.(f)) .* S2_81_hh_indicator.CVA_vars.irf.(f); 
    S2_81_hh_indicator.aLP_RS_CVA.mu.(f) = mean(S2_81_hh_indicator.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.aLP_RS_CVA.band_lo.(f) = prctile(S2_81_hh_indicator.aLP_RS_CVA.irf.(f),pct_lo,2);
    S2_81_hh_indicator.aLP_RS_CVA.band_hi.(f) = prctile(S2_81_hh_indicator.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S2_81_hh_indicator.aLP_MSPE_GMA.irf.(f) = S2_81_hh_indicator.aLP_MSPE.alpha.(f) .* S2_81_hh_indicator.GMA_lps.irf.(f)...
    + (1-S2_81_hh_indicator.aLP_MSPE.alpha.(f)) .* S2_81_hh_indicator.GMA_vars.irf.(f); 
    S2_81_hh_indicator.aLP_MSPE_GMA.mu.(f) = mean(S2_81_hh_indicator.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.aLP_MSPE_GMA.band_lo.(f) = prctile(S2_81_hh_indicator.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S2_81_hh_indicator.aLP_MSPE_GMA.band_hi.(f) = prctile(S2_81_hh_indicator.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S2_81_hh_indicator.aLP_MSPE_CVA.irf.(f) = S2_81_hh_indicator.aLP_MSPE.alpha.(f) .* S2_81_hh_indicator.CVA_lps.irf.(f)...
    + (1-S2_81_hh_indicator.aLP_MSPE.alpha.(f)) .* S2_81_hh_indicator.CVA_vars.irf.(f); 
    S2_81_hh_indicator.aLP_MSPE_CVA.mu.(f) = mean(S2_81_hh_indicator.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.aLP_MSPE_CVA.band_lo.(f) = prctile(S2_81_hh_indicator.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S2_81_hh_indicator.aLP_MSPE_CVA.band_hi.(f) = prctile(S2_81_hh_indicator.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S2_81_hh_indicator.aLP_RS_GMA_girf.irf.(f) = S2_81_hh_indicator.aLP_RS.alpha.(f) .* S2_81_hh_indicator.GMA_lps.irf.(f)...
    + (1-S2_81_hh_indicator.aLP_RS.alpha.(f)) .* S2_81_hh_indicator.GMA_vars_girf.irf.(f); 
    S2_81_hh_indicator.aLP_RS_GMA_girf.mu.(f) = mean(S2_81_hh_indicator.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.aLP_RS_GMA_girf.band_lo.(f) = prctile(S2_81_hh_indicator.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S2_81_hh_indicator.aLP_RS_GMA_girf.band_hi.(f) = prctile(S2_81_hh_indicator.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_81_hh_indicator.aLP_RS_CVA_girf.irf.(f) = S2_81_hh_indicator.aLP_RS.alpha.(f) .* S2_81_hh_indicator.CVA_lps.irf.(f)...
    + (1-S2_81_hh_indicator.aLP_RS.alpha.(f)) .* S2_81_hh_indicator.CVA_vars_girf.irf.(f); 
    S2_81_hh_indicator.aLP_RS_CVA_girf.mu.(f) = mean(S2_81_hh_indicator.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.aLP_RS_CVA_girf.band_lo.(f) = prctile(S2_81_hh_indicator.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S2_81_hh_indicator.aLP_RS_CVA_girf.band_hi.(f) = prctile(S2_81_hh_indicator.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S2_81_hh_indicator.aLP_MSPE_GMA_girf.irf.(f) = S2_81_hh_indicator.aLP_MSPE.alpha.(f) .* S2_81_hh_indicator.GMA_lps.irf.(f)...
    + (1-S2_81_hh_indicator.aLP_MSPE.alpha.(f)) .* S2_81_hh_indicator.GMA_vars_girf.irf.(f); 
    S2_81_hh_indicator.aLP_MSPE_GMA_girf.mu.(f) = mean(S2_81_hh_indicator.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S2_81_hh_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S2_81_hh_indicator.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S2_81_hh_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_81_hh_indicator.aLP_MSPE_CVA_girf.irf.(f) = S2_81_hh_indicator.aLP_MSPE.alpha.(f) .* S2_81_hh_indicator.CVA_lps.irf.(f)...
    + (1-S2_81_hh_indicator.aLP_MSPE.alpha.(f)) .* S2_81_hh_indicator.CVA_vars_girf.irf.(f); 
    S2_81_hh_indicator.aLP_MSPE_CVA_girf.mu.(f) = mean(S2_81_hh_indicator.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S2_81_hh_indicator.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S2_81_hh_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S2_81_hh_indicator.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S2_81_hh_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S2_81_hh_indicator','-append')
%% ------------------------------------------------------
% 2-3. SPF IQR from 1981
% S2_81_pro_indcator
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\indicator\lag12\Bloom_S2_81_pro_indicator.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S2_81_pro_indicator.(names{i}).irf.(f) = irfpaths.(names{i}).S2_81_pro.(f);
        S2_81_pro_indicator.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S2_81_pro.(f),2);
        S2_81_pro_indicator.(names{i}).band_hi.(f) = band_hi.(names{i}).S2_81_pro.(f);
        S2_81_pro_indicator.(names{i}).band_lo.(f) = band_lo.(names{i}).S2_81_pro.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S2_81_pro.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S2_81_pro_indicator.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_81_pro.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S2_81_pro_indicator.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S2_81_pro.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S2_81_pro_indicator.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S2_81_pro.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S2_81_pro_indicator.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S2_81_pro.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S2_81_pro_indicator.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S2_81_pro.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S2_81_pro_indicator.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_indicator.GMA_lps.irf.(f) = wirf_val;% H x B
    S2_81_pro_indicator.GMA_lps.mu.(f) = mean(S2_81_pro_indicator.GMA_lps.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_indicator.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S2_81_pro_indicator.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_indicator.CVA_lps.irf.(f) = wirf_val;% H x B
    S2_81_pro_indicator.CVA_lps.mu.(f) = mean(S2_81_pro_indicator.CVA_lps.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_indicator.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S2_81_pro.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_81_pro_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_indicator.GMA_vars.irf.(f) = wirf_val;% H x B
    S2_81_pro_indicator.GMA_vars.mu.(f) = mean(S2_81_pro_indicator.GMA_vars.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_indicator.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_81_pro_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_indicator.CVA_vars.irf.(f) = wirf_val;
    S2_81_pro_indicator.CVA_vars.mu.(f) = mean(S2_81_pro_indicator.CVA_vars.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_indicator.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S2_81_pro.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S2_81_pro_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_indicator.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S2_81_pro_indicator.GMA_vars_girf.mu.(f) = mean(S2_81_pro_indicator.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_indicator.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S2_81_pro_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S2_81_pro_indicator.CVA_vars_girf.irf.(f) = wirf_val;
    S2_81_pro_indicator.CVA_vars_girf.mu.(f) = mean(S2_81_pro_indicator.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S2_81_pro_indicator.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S2_81_pro.state0_vol_empl = tss.S2_81_pro.state0_lps_empl;
    tss.S2_81_pro.state0_ff_empl = tss.S2_81_pro.state0_lps_empl;
    tss.S2_81_pro.state1_vol_empl = tss.S2_81_pro.state1_lps_empl;
    tss.S2_81_pro.state1_ff_empl = tss.S2_81_pro.state1_lps_empl;
    tss.S2_81_pro.state0_vol_ip = tss.S2_81_pro.state0_lps_ip;
    tss.S2_81_pro.state0_ff_ip = tss.S2_81_pro.state0_lps_ip;
    tss.S2_81_pro.state1_vol_ip = tss.S2_81_pro.state1_lps_ip;
    tss.S2_81_pro.state1_ff_ip = tss.S2_81_pro.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S2_81_pro.(f) ./ tss.S2_81_pro.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S2_81_pro.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S2_81_pro.state0_vol_empl = tss.S2_81_pro.state0_vars_empl;
    tss.S2_81_pro.state0_ff_empl = tss.S2_81_pro.state0_vars_empl;
    tss.S2_81_pro.state1_vol_empl = tss.S2_81_pro.state1_vars_empl;
    tss.S2_81_pro.state1_ff_empl = tss.S2_81_pro.state1_vars_empl;
    tss.S2_81_pro.state0_vol_ip = tss.S2_81_pro.state0_vars_ip;
    tss.S2_81_pro.state0_ff_ip = tss.S2_81_pro.state0_vars_ip;
    tss.S2_81_pro.state1_vol_ip = tss.S2_81_pro.state1_vars_ip;
    tss.S2_81_pro.state1_ff_ip = tss.S2_81_pro.state1_vars_ip;
    tss.S2_81_pro.linear_vol_empl = tss.S2_81_pro.linear_vars_empl;
    tss.S2_81_pro.linear_vol_ip = tss.S2_81_pro.linear_vars_ip;
    tss.S2_81_pro.linear_ff_empl = tss.S2_81_pro.linear_vars_empl;
    tss.S2_81_pro.linear_ff_ip = tss.S2_81_pro.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S2_81_pro.(f) ./ tss.S2_81_pro.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S2_81_pro.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S2_81_pro_indicator.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S2_81_pro_indicator.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S2_81_pro_indicator.aLP_RS.alpha.(f) .* S2_81_pro_indicator.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S2_81_pro_indicator.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_81_pro_indicator.GMA_vars.weight.(f)'; % H x 2
    S2_81_pro_indicator.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S2_81_pro_indicator.aLP_RS.alpha.(f) .* S2_81_pro_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_81_pro_indicator.CVA_vars.weight.(f)';
    S2_81_pro_indicator.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S2_81_pro_indicator.aLP_MSPE.alpha.(f) .* S2_81_pro_indicator.GMA_lps.weight.(f)';
    alpha_var = 1 - S2_81_pro_indicator.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S2_81_pro_indicator.GMA_vars.weight.(f)';
    S2_81_pro_indicator.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S2_81_pro_indicator.aLP_MSPE.alpha.(f) .* S2_81_pro_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S2_81_pro_indicator.CVA_vars.weight.(f)';
    S2_81_pro_indicator.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S2_81_pro_indicator.aLP_RS_GMA.irf.(f) = S2_81_pro_indicator.aLP_RS.alpha.(f) .* S2_81_pro_indicator.GMA_lps.irf.(f)...
    + (1-S2_81_pro_indicator.aLP_RS.alpha.(f)) .* S2_81_pro_indicator.GMA_vars.irf.(f); 
    S2_81_pro_indicator.aLP_RS_GMA.mu.(f) = mean(S2_81_pro_indicator.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.aLP_RS_GMA.band_lo.(f) = prctile(S2_81_pro_indicator.aLP_RS_GMA.irf.(f),pct_lo,2);
    S2_81_pro_indicator.aLP_RS_GMA.band_hi.(f) = prctile(S2_81_pro_indicator.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S2_81_pro_indicator.aLP_RS_CVA.irf.(f) = S2_81_pro_indicator.aLP_RS.alpha.(f) .* S2_81_pro_indicator.CVA_lps.irf.(f)...
    + (1-S2_81_pro_indicator.aLP_RS.alpha.(f)) .* S2_81_pro_indicator.CVA_vars.irf.(f); 
    S2_81_pro_indicator.aLP_RS_CVA.mu.(f) = mean(S2_81_pro_indicator.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.aLP_RS_CVA.band_lo.(f) = prctile(S2_81_pro_indicator.aLP_RS_CVA.irf.(f),pct_lo,2);
    S2_81_pro_indicator.aLP_RS_CVA.band_hi.(f) = prctile(S2_81_pro_indicator.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S2_81_pro_indicator.aLP_MSPE_GMA.irf.(f) = S2_81_pro_indicator.aLP_MSPE.alpha.(f) .* S2_81_pro_indicator.GMA_lps.irf.(f)...
    + (1-S2_81_pro_indicator.aLP_MSPE.alpha.(f)) .* S2_81_pro_indicator.GMA_vars.irf.(f); 
    S2_81_pro_indicator.aLP_MSPE_GMA.mu.(f) = mean(S2_81_pro_indicator.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.aLP_MSPE_GMA.band_lo.(f) = prctile(S2_81_pro_indicator.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S2_81_pro_indicator.aLP_MSPE_GMA.band_hi.(f) = prctile(S2_81_pro_indicator.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S2_81_pro_indicator.aLP_MSPE_CVA.irf.(f) = S2_81_pro_indicator.aLP_MSPE.alpha.(f) .* S2_81_pro_indicator.CVA_lps.irf.(f)...
    + (1-S2_81_pro_indicator.aLP_MSPE.alpha.(f)) .* S2_81_pro_indicator.CVA_vars.irf.(f); 
    S2_81_pro_indicator.aLP_MSPE_CVA.mu.(f) = mean(S2_81_pro_indicator.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.aLP_MSPE_CVA.band_lo.(f) = prctile(S2_81_pro_indicator.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S2_81_pro_indicator.aLP_MSPE_CVA.band_hi.(f) = prctile(S2_81_pro_indicator.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S2_81_pro_indicator.aLP_RS_GMA_girf.irf.(f) = S2_81_pro_indicator.aLP_RS.alpha.(f) .* S2_81_pro_indicator.GMA_lps.irf.(f)...
    + (1-S2_81_pro_indicator.aLP_RS.alpha.(f)) .* S2_81_pro_indicator.GMA_vars_girf.irf.(f); 
    S2_81_pro_indicator.aLP_RS_GMA_girf.mu.(f) = mean(S2_81_pro_indicator.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.aLP_RS_GMA_girf.band_lo.(f) = prctile(S2_81_pro_indicator.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S2_81_pro_indicator.aLP_RS_GMA_girf.band_hi.(f) = prctile(S2_81_pro_indicator.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_81_pro_indicator.aLP_RS_CVA_girf.irf.(f) = S2_81_pro_indicator.aLP_RS.alpha.(f) .* S2_81_pro_indicator.CVA_lps.irf.(f)...
    + (1-S2_81_pro_indicator.aLP_RS.alpha.(f)) .* S2_81_pro_indicator.CVA_vars_girf.irf.(f); 
    S2_81_pro_indicator.aLP_RS_CVA_girf.mu.(f) = mean(S2_81_pro_indicator.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.aLP_RS_CVA_girf.band_lo.(f) = prctile(S2_81_pro_indicator.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S2_81_pro_indicator.aLP_RS_CVA_girf.band_hi.(f) = prctile(S2_81_pro_indicator.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S2_81_pro_indicator.aLP_MSPE_GMA_girf.irf.(f) = S2_81_pro_indicator.aLP_MSPE.alpha.(f) .* S2_81_pro_indicator.GMA_lps.irf.(f)...
    + (1-S2_81_pro_indicator.aLP_MSPE.alpha.(f)) .* S2_81_pro_indicator.GMA_vars_girf.irf.(f); 
    S2_81_pro_indicator.aLP_MSPE_GMA_girf.mu.(f) = mean(S2_81_pro_indicator.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S2_81_pro_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S2_81_pro_indicator.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S2_81_pro_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S2_81_pro_indicator.aLP_MSPE_CVA_girf.irf.(f) = S2_81_pro_indicator.aLP_MSPE.alpha.(f) .* S2_81_pro_indicator.CVA_lps.irf.(f)...
    + (1-S2_81_pro_indicator.aLP_MSPE.alpha.(f)) .* S2_81_pro_indicator.CVA_vars_girf.irf.(f); 
    S2_81_pro_indicator.aLP_MSPE_CVA_girf.mu.(f) = mean(S2_81_pro_indicator.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S2_81_pro_indicator.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S2_81_pro_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S2_81_pro_indicator.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S2_81_pro_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S2_81_pro_indicator','-append')
%% -------------------------------------------------------
% 3. Low vs. High inflation volatility
% 3-1. Fed target vs actual inflation rate by abolute
% S3_abs_indcator
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\indicator\lag12\Bloom_S3_abs_indicator.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S3_abs_indicator.(names{i}).irf.(f) = irfpaths.(names{i}).S3_abs.(f);
        S3_abs_indicator.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S3_abs.(f),2);
        S3_abs_indicator.(names{i}).band_hi.(f) = band_hi.(names{i}).S3_abs.(f);
        S3_abs_indicator.(names{i}).band_lo.(f) = band_lo.(names{i}).S3_abs.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S3_abs.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S3_abs_indicator.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_abs.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S3_abs_indicator.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_abs.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S3_abs_indicator.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S3_abs.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S3_abs_indicator.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S3_abs.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S3_abs_indicator.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S3_abs.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S3_abs_indicator.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_indicator.GMA_lps.irf.(f) = wirf_val;% H x B
    S3_abs_indicator.GMA_lps.mu.(f) = mean(S3_abs_indicator.GMA_lps.irf.(f), 2, "omitmissing");
    S3_abs_indicator.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_indicator.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S3_abs_indicator.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_indicator.CVA_lps.irf.(f) = wirf_val;% H x B
    S3_abs_indicator.CVA_lps.mu.(f) = mean(S3_abs_indicator.CVA_lps.irf.(f), 2, "omitmissing");
    S3_abs_indicator.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_indicator.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S3_abs.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_abs_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_indicator.GMA_vars.irf.(f) = wirf_val;% H x B
    S3_abs_indicator.GMA_vars.mu.(f) = mean(S3_abs_indicator.GMA_vars.irf.(f), 2, "omitmissing");
    S3_abs_indicator.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_indicator.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_abs_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_indicator.CVA_vars.irf.(f) = wirf_val;
    S3_abs_indicator.CVA_vars.mu.(f) = mean(S3_abs_indicator.CVA_vars.irf.(f), 2, "omitmissing");
    S3_abs_indicator.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_indicator.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S3_abs.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_abs_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_indicator.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S3_abs_indicator.GMA_vars_girf.mu.(f) = mean(S3_abs_indicator.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S3_abs_indicator.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_indicator.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_abs_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_abs_indicator.CVA_vars_girf.irf.(f) = wirf_val;
    S3_abs_indicator.CVA_vars_girf.mu.(f) = mean(S3_abs_indicator.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S3_abs_indicator.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_abs_indicator.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S3_abs.state0_vol_empl = tss.S3_abs.state0_lps_empl;
    tss.S3_abs.state0_ff_empl = tss.S3_abs.state0_lps_empl;
    tss.S3_abs.state1_vol_empl = tss.S3_abs.state1_lps_empl;
    tss.S3_abs.state1_ff_empl = tss.S3_abs.state1_lps_empl;
    tss.S3_abs.state0_vol_ip = tss.S3_abs.state0_lps_ip;
    tss.S3_abs.state0_ff_ip = tss.S3_abs.state0_lps_ip;
    tss.S3_abs.state1_vol_ip = tss.S3_abs.state1_lps_ip;
    tss.S3_abs.state1_ff_ip = tss.S3_abs.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S3_abs.(f) ./ tss.S3_abs.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S3_abs.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S3_abs.state0_vol_empl = tss.S3_abs.state0_vars_empl;
    tss.S3_abs.state0_ff_empl = tss.S3_abs.state0_vars_empl;
    tss.S3_abs.state1_vol_empl = tss.S3_abs.state1_vars_empl;
    tss.S3_abs.state1_ff_empl = tss.S3_abs.state1_vars_empl;
    tss.S3_abs.state0_vol_ip = tss.S3_abs.state0_vars_ip;
    tss.S3_abs.state0_ff_ip = tss.S3_abs.state0_vars_ip;
    tss.S3_abs.state1_vol_ip = tss.S3_abs.state1_vars_ip;
    tss.S3_abs.state1_ff_ip = tss.S3_abs.state1_vars_ip;
    tss.S3_abs.linear_vol_empl = tss.S3_abs.linear_vars_empl;
    tss.S3_abs.linear_vol_ip = tss.S3_abs.linear_vars_ip;
    tss.S3_abs.linear_ff_empl = tss.S3_abs.linear_vars_empl;
    tss.S3_abs.linear_ff_ip = tss.S3_abs.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S3_abs.(f) ./ tss.S3_abs.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S3_abs.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S3_abs_indicator.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S3_abs_indicator.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S3_abs_indicator.aLP_RS.alpha.(f) .* S3_abs_indicator.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S3_abs_indicator.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_abs_indicator.GMA_vars.weight.(f)'; % H x 2
    S3_abs_indicator.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S3_abs_indicator.aLP_RS.alpha.(f) .* S3_abs_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_abs_indicator.CVA_vars.weight.(f)';
    S3_abs_indicator.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S3_abs_indicator.aLP_MSPE.alpha.(f) .* S3_abs_indicator.GMA_lps.weight.(f)';
    alpha_var = 1 - S3_abs_indicator.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_abs_indicator.GMA_vars.weight.(f)';
    S3_abs_indicator.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S3_abs_indicator.aLP_MSPE.alpha.(f) .* S3_abs_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_abs_indicator.CVA_vars.weight.(f)';
    S3_abs_indicator.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S3_abs_indicator.aLP_RS_GMA.irf.(f) = S3_abs_indicator.aLP_RS.alpha.(f) .* S3_abs_indicator.GMA_lps.irf.(f)...
    + (1-S3_abs_indicator.aLP_RS.alpha.(f)) .* S3_abs_indicator.GMA_vars.irf.(f); 
    S3_abs_indicator.aLP_RS_GMA.mu.(f) = mean(S3_abs_indicator.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S3_abs_indicator.aLP_RS_GMA.band_lo.(f) = prctile(S3_abs_indicator.aLP_RS_GMA.irf.(f),pct_lo,2);
    S3_abs_indicator.aLP_RS_GMA.band_hi.(f) = prctile(S3_abs_indicator.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S3_abs_indicator.aLP_RS_CVA.irf.(f) = S3_abs_indicator.aLP_RS.alpha.(f) .* S3_abs_indicator.CVA_lps.irf.(f)...
    + (1-S3_abs_indicator.aLP_RS.alpha.(f)) .* S3_abs_indicator.CVA_vars.irf.(f); 
    S3_abs_indicator.aLP_RS_CVA.mu.(f) = mean(S3_abs_indicator.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S3_abs_indicator.aLP_RS_CVA.band_lo.(f) = prctile(S3_abs_indicator.aLP_RS_CVA.irf.(f),pct_lo,2);
    S3_abs_indicator.aLP_RS_CVA.band_hi.(f) = prctile(S3_abs_indicator.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S3_abs_indicator.aLP_MSPE_GMA.irf.(f) = S3_abs_indicator.aLP_MSPE.alpha.(f) .* S3_abs_indicator.GMA_lps.irf.(f)...
    + (1-S3_abs_indicator.aLP_MSPE.alpha.(f)) .* S3_abs_indicator.GMA_vars.irf.(f); 
    S3_abs_indicator.aLP_MSPE_GMA.mu.(f) = mean(S3_abs_indicator.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S3_abs_indicator.aLP_MSPE_GMA.band_lo.(f) = prctile(S3_abs_indicator.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S3_abs_indicator.aLP_MSPE_GMA.band_hi.(f) = prctile(S3_abs_indicator.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S3_abs_indicator.aLP_MSPE_CVA.irf.(f) = S3_abs_indicator.aLP_MSPE.alpha.(f) .* S3_abs_indicator.CVA_lps.irf.(f)...
    + (1-S3_abs_indicator.aLP_MSPE.alpha.(f)) .* S3_abs_indicator.CVA_vars.irf.(f); 
    S3_abs_indicator.aLP_MSPE_CVA.mu.(f) = mean(S3_abs_indicator.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S3_abs_indicator.aLP_MSPE_CVA.band_lo.(f) = prctile(S3_abs_indicator.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S3_abs_indicator.aLP_MSPE_CVA.band_hi.(f) = prctile(S3_abs_indicator.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S3_abs_indicator.aLP_RS_GMA_girf.irf.(f) = S3_abs_indicator.aLP_RS.alpha.(f) .* S3_abs_indicator.GMA_lps.irf.(f)...
    + (1-S3_abs_indicator.aLP_RS.alpha.(f)) .* S3_abs_indicator.GMA_vars_girf.irf.(f); 
    S3_abs_indicator.aLP_RS_GMA_girf.mu.(f) = mean(S3_abs_indicator.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S3_abs_indicator.aLP_RS_GMA_girf.band_lo.(f) = prctile(S3_abs_indicator.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S3_abs_indicator.aLP_RS_GMA_girf.band_hi.(f) = prctile(S3_abs_indicator.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_abs_indicator.aLP_RS_CVA_girf.irf.(f) = S3_abs_indicator.aLP_RS.alpha.(f) .* S3_abs_indicator.CVA_lps.irf.(f)...
    + (1-S3_abs_indicator.aLP_RS.alpha.(f)) .* S3_abs_indicator.CVA_vars_girf.irf.(f); 
    S3_abs_indicator.aLP_RS_CVA_girf.mu.(f) = mean(S3_abs_indicator.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S3_abs_indicator.aLP_RS_CVA_girf.band_lo.(f) = prctile(S3_abs_indicator.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S3_abs_indicator.aLP_RS_CVA_girf.band_hi.(f) = prctile(S3_abs_indicator.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S3_abs_indicator.aLP_MSPE_GMA_girf.irf.(f) = S3_abs_indicator.aLP_MSPE.alpha.(f) .* S3_abs_indicator.GMA_lps.irf.(f)...
    + (1-S3_abs_indicator.aLP_MSPE.alpha.(f)) .* S3_abs_indicator.GMA_vars_girf.irf.(f); 
    S3_abs_indicator.aLP_MSPE_GMA_girf.mu.(f) = mean(S3_abs_indicator.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S3_abs_indicator.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S3_abs_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S3_abs_indicator.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S3_abs_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_abs_indicator.aLP_MSPE_CVA_girf.irf.(f) = S3_abs_indicator.aLP_MSPE.alpha.(f) .* S3_abs_indicator.CVA_lps.irf.(f)...
    + (1-S3_abs_indicator.aLP_MSPE.alpha.(f)) .* S3_abs_indicator.CVA_vars_girf.irf.(f); 
    S3_abs_indicator.aLP_MSPE_CVA_girf.mu.(f) = mean(S3_abs_indicator.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S3_abs_indicator.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S3_abs_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S3_abs_indicator.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S3_abs_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S3_abs_indicator','-append')
%% ------------------------------------------------------
% 3-2. Fed target vs actual inflation rate by 75th percentile
% S3_rel_indcator
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\indicator\lag12\Bloom_S3_rel_indicator.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S3_rel_indicator.(names{i}).irf.(f) = irfpaths.(names{i}).S3_rel.(f);
        S3_rel_indicator.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S3_rel.(f),2);
        S3_rel_indicator.(names{i}).band_hi.(f) = band_hi.(names{i}).S3_rel.(f);
        S3_rel_indicator.(names{i}).band_lo.(f) = band_lo.(names{i}).S3_rel.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S3_rel.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S3_rel_indicator.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_rel.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S3_rel_indicator.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_rel.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S3_rel_indicator.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S3_rel.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S3_rel_indicator.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S3_rel.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S3_rel_indicator.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S3_rel.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S3_rel_indicator.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_indicator.GMA_lps.irf.(f) = wirf_val;% H x B
    S3_rel_indicator.GMA_lps.mu.(f) = mean(S3_rel_indicator.GMA_lps.irf.(f), 2, "omitmissing");
    S3_rel_indicator.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_indicator.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S3_rel_indicator.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_indicator.CVA_lps.irf.(f) = wirf_val;% H x B
    S3_rel_indicator.CVA_lps.mu.(f) = mean(S3_rel_indicator.CVA_lps.irf.(f), 2, "omitmissing");
    S3_rel_indicator.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_indicator.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S3_rel.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_rel_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_indicator.GMA_vars.irf.(f) = wirf_val;% H x B
    S3_rel_indicator.GMA_vars.mu.(f) = mean(S3_rel_indicator.GMA_vars.irf.(f), 2, "omitmissing");
    S3_rel_indicator.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_indicator.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_rel_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_indicator.CVA_vars.irf.(f) = wirf_val;
    S3_rel_indicator.CVA_vars.mu.(f) = mean(S3_rel_indicator.CVA_vars.irf.(f), 2, "omitmissing");
    S3_rel_indicator.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_indicator.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S3_rel.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_rel_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_indicator.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S3_rel_indicator.GMA_vars_girf.mu.(f) = mean(S3_rel_indicator.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S3_rel_indicator.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_indicator.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_rel_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_rel_indicator.CVA_vars_girf.irf.(f) = wirf_val;
    S3_rel_indicator.CVA_vars_girf.mu.(f) = mean(S3_rel_indicator.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S3_rel_indicator.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_rel_indicator.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S3_rel.state0_vol_empl = tss.S3_rel.state0_lps_empl;
    tss.S3_rel.state0_ff_empl = tss.S3_rel.state0_lps_empl;
    tss.S3_rel.state1_vol_empl = tss.S3_rel.state1_lps_empl;
    tss.S3_rel.state1_ff_empl = tss.S3_rel.state1_lps_empl;
    tss.S3_rel.state0_vol_ip = tss.S3_rel.state0_lps_ip;
    tss.S3_rel.state0_ff_ip = tss.S3_rel.state0_lps_ip;
    tss.S3_rel.state1_vol_ip = tss.S3_rel.state1_lps_ip;
    tss.S3_rel.state1_ff_ip = tss.S3_rel.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S3_rel.(f) ./ tss.S3_rel.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S3_rel.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S3_rel.state0_vol_empl = tss.S3_rel.state0_vars_empl;
    tss.S3_rel.state0_ff_empl = tss.S3_rel.state0_vars_empl;
    tss.S3_rel.state1_vol_empl = tss.S3_rel.state1_vars_empl;
    tss.S3_rel.state1_ff_empl = tss.S3_rel.state1_vars_empl;
    tss.S3_rel.state0_vol_ip = tss.S3_rel.state0_vars_ip;
    tss.S3_rel.state0_ff_ip = tss.S3_rel.state0_vars_ip;
    tss.S3_rel.state1_vol_ip = tss.S3_rel.state1_vars_ip;
    tss.S3_rel.state1_ff_ip = tss.S3_rel.state1_vars_ip;
    tss.S3_rel.linear_vol_empl = tss.S3_rel.linear_vars_empl;
    tss.S3_rel.linear_vol_ip = tss.S3_rel.linear_vars_ip;
    tss.S3_rel.linear_ff_empl = tss.S3_rel.linear_vars_empl;
    tss.S3_rel.linear_ff_ip = tss.S3_rel.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S3_rel.(f) ./ tss.S3_rel.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S3_rel.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S3_rel_indicator.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S3_rel_indicator.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S3_rel_indicator.aLP_RS.alpha.(f) .* S3_rel_indicator.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S3_rel_indicator.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_rel_indicator.GMA_vars.weight.(f)'; % H x 2
    S3_rel_indicator.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S3_rel_indicator.aLP_RS.alpha.(f) .* S3_rel_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_rel_indicator.CVA_vars.weight.(f)';
    S3_rel_indicator.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S3_rel_indicator.aLP_MSPE.alpha.(f) .* S3_rel_indicator.GMA_lps.weight.(f)';
    alpha_var = 1 - S3_rel_indicator.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_rel_indicator.GMA_vars.weight.(f)';
    S3_rel_indicator.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S3_rel_indicator.aLP_MSPE.alpha.(f) .* S3_rel_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_rel_indicator.CVA_vars.weight.(f)';
    S3_rel_indicator.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S3_rel_indicator.aLP_RS_GMA.irf.(f) = S3_rel_indicator.aLP_RS.alpha.(f) .* S3_rel_indicator.GMA_lps.irf.(f)...
    + (1-S3_rel_indicator.aLP_RS.alpha.(f)) .* S3_rel_indicator.GMA_vars.irf.(f); 
    S3_rel_indicator.aLP_RS_GMA.mu.(f) = mean(S3_rel_indicator.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S3_rel_indicator.aLP_RS_GMA.band_lo.(f) = prctile(S3_rel_indicator.aLP_RS_GMA.irf.(f),pct_lo,2);
    S3_rel_indicator.aLP_RS_GMA.band_hi.(f) = prctile(S3_rel_indicator.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S3_rel_indicator.aLP_RS_CVA.irf.(f) = S3_rel_indicator.aLP_RS.alpha.(f) .* S3_rel_indicator.CVA_lps.irf.(f)...
    + (1-S3_rel_indicator.aLP_RS.alpha.(f)) .* S3_rel_indicator.CVA_vars.irf.(f); 
    S3_rel_indicator.aLP_RS_CVA.mu.(f) = mean(S3_rel_indicator.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S3_rel_indicator.aLP_RS_CVA.band_lo.(f) = prctile(S3_rel_indicator.aLP_RS_CVA.irf.(f),pct_lo,2);
    S3_rel_indicator.aLP_RS_CVA.band_hi.(f) = prctile(S3_rel_indicator.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S3_rel_indicator.aLP_MSPE_GMA.irf.(f) = S3_rel_indicator.aLP_MSPE.alpha.(f) .* S3_rel_indicator.GMA_lps.irf.(f)...
    + (1-S3_rel_indicator.aLP_MSPE.alpha.(f)) .* S3_rel_indicator.GMA_vars.irf.(f); 
    S3_rel_indicator.aLP_MSPE_GMA.mu.(f) = mean(S3_rel_indicator.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S3_rel_indicator.aLP_MSPE_GMA.band_lo.(f) = prctile(S3_rel_indicator.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S3_rel_indicator.aLP_MSPE_GMA.band_hi.(f) = prctile(S3_rel_indicator.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S3_rel_indicator.aLP_MSPE_CVA.irf.(f) = S3_rel_indicator.aLP_MSPE.alpha.(f) .* S3_rel_indicator.CVA_lps.irf.(f)...
    + (1-S3_rel_indicator.aLP_MSPE.alpha.(f)) .* S3_rel_indicator.CVA_vars.irf.(f); 
    S3_rel_indicator.aLP_MSPE_CVA.mu.(f) = mean(S3_rel_indicator.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S3_rel_indicator.aLP_MSPE_CVA.band_lo.(f) = prctile(S3_rel_indicator.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S3_rel_indicator.aLP_MSPE_CVA.band_hi.(f) = prctile(S3_rel_indicator.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S3_rel_indicator.aLP_RS_GMA_girf.irf.(f) = S3_rel_indicator.aLP_RS.alpha.(f) .* S3_rel_indicator.GMA_lps.irf.(f)...
    + (1-S3_rel_indicator.aLP_RS.alpha.(f)) .* S3_rel_indicator.GMA_vars_girf.irf.(f); 
    S3_rel_indicator.aLP_RS_GMA_girf.mu.(f) = mean(S3_rel_indicator.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S3_rel_indicator.aLP_RS_GMA_girf.band_lo.(f) = prctile(S3_rel_indicator.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S3_rel_indicator.aLP_RS_GMA_girf.band_hi.(f) = prctile(S3_rel_indicator.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_rel_indicator.aLP_RS_CVA_girf.irf.(f) = S3_rel_indicator.aLP_RS.alpha.(f) .* S3_rel_indicator.CVA_lps.irf.(f)...
    + (1-S3_rel_indicator.aLP_RS.alpha.(f)) .* S3_rel_indicator.CVA_vars_girf.irf.(f); 
    S3_rel_indicator.aLP_RS_CVA_girf.mu.(f) = mean(S3_rel_indicator.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S3_rel_indicator.aLP_RS_CVA_girf.band_lo.(f) = prctile(S3_rel_indicator.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S3_rel_indicator.aLP_RS_CVA_girf.band_hi.(f) = prctile(S3_rel_indicator.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S3_rel_indicator.aLP_MSPE_GMA_girf.irf.(f) = S3_rel_indicator.aLP_MSPE.alpha.(f) .* S3_rel_indicator.GMA_lps.irf.(f)...
    + (1-S3_rel_indicator.aLP_MSPE.alpha.(f)) .* S3_rel_indicator.GMA_vars_girf.irf.(f); 
    S3_rel_indicator.aLP_MSPE_GMA_girf.mu.(f) = mean(S3_rel_indicator.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S3_rel_indicator.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S3_rel_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S3_rel_indicator.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S3_rel_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_rel_indicator.aLP_MSPE_CVA_girf.irf.(f) = S3_rel_indicator.aLP_MSPE.alpha.(f) .* S3_rel_indicator.CVA_lps.irf.(f)...
    + (1-S3_rel_indicator.aLP_MSPE.alpha.(f)) .* S3_rel_indicator.CVA_vars_girf.irf.(f); 
    S3_rel_indicator.aLP_MSPE_CVA_girf.mu.(f) = mean(S3_rel_indicator.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S3_rel_indicator.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S3_rel_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S3_rel_indicator.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S3_rel_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S3_rel_indicator','-append')
%% ------------------------------------------------------
% 3-3. by rolling sigma 75th percentile
% S3_75thr_indcator
% -------------------------------------------------------
clc
clear all 
load('Combine_Results_Emp.mat')
load("Results_Empirical\indicator\lag12\Bloom_S3_75thr_indicator.mat")
% settings 
names  = {'svar', 'bvar', 'lp', 'lp_contmp', 'girf'};
fields = {'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip', 'state1_ff_ip'}; 
% sigle estimators: irf, bands (68% ci)
for i = 1:length(names)
    for j = 1:length(fields)
        f = fields{j};
        S3_75thr_indicator.(names{i}).irf.(f) = irfpaths.(names{i}).S3_75thr.(f);
        S3_75thr_indicator.(names{i}).mu.(f) = mean(irfpaths.(names{i}).S3_75thr.(f),2);
        S3_75thr_indicator.(names{i}).band_hi.(f) = band_hi.(names{i}).S3_75thr.(f);
        S3_75thr_indicator.(names{i}).band_lo.(f) = band_lo.(names{i}).S3_75thr.(f);
    end
end
% --------------------------------------------------------------------
% mavg: gma, cva 
% --------------------------------------------------------------------
% GMA ----------------------------------------------------
%{
% comb_fields ={'state0_vol_empl', 'state1_vol_empl', 'state0_ff_empl', 'state1_ff_empl',...
    'state0_vol_ip', 'state1_vol_ip', 'state0_ff_ip','state1_ff_ip',...
    'linear_vol_empl','linear_ff_empl','linear_vol_ip','linear_ff_ip'};
%}
lp_names  = {'lp', 'lp_contmp'};
var_names = {'svar', 'bvar'};
var_names_rob ={'girf','bvar'};
% lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(lp_names{i}).S3_75thr.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % n_models x H 
    S3_75thr_indicator.GMA_lps.weight.(f) = gma_val; % n_models x H 
end
% vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 

    % Collect raw BICs first into n_models x H matrix
    GIC_mat = zeros(n_models, settings.est.IRF_hor);
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_75thr.(f); % H x 1
        GIC_mat(i,:) = GIC_val(:)';                 % force row
    end
    % Subtract column minima (log-sum-exp trick) to prevent NaN
    minBIC = min(GIC_mat,[],1);                       % 1 x H
    weights_raw = exp(-(GIC_mat - minBIC) / 2);       % n_models x H
    % Normalize across models
    gma_val = weights_raw ./ sum(weights_raw,1);      % n_models x H
    S3_75thr_indicator.GMA_vars.weight.(f) = gma_val; 
end
%}
%{
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j}; 
    GIC_row = zeros(n_models, settings.est.IRF_hor); % [2 x H ]
    for i = 1:n_models
        GIC_val = gic.(var_names{i}).S3_75thr.(f); % H x 1
        GIC_row(i,:) = exp(-GIC_val/2); % 
    end
    % sum/normalize for this field only
    GIC_row_sum = sum(GIC_row,1);        % 1 x H 
    gma_val     = GIC_row ./ GIC_row_sum; % 2 x H 
    S3_75thr_indicator.GMA_vars.weight.(f) = gma_val; % 2 x H 
end
%}
%
clear GIC_row GIC_row_sum GIC_val gma_val GIC_mat minBIC weights_raw 
% vars-girf share the same weight with vars
% CVA -----------------------------------------------
for j = 1:length(fields)
    f = fields{j};
    %
    cva_val = CVA_w_lps.S3_75thr.(f); % H x 2
    cva_val_trans = cva_val'; % 2 x H 
    S3_75thr_indicator.CVA_lps.weight.(f) = cva_val_trans;
    %
    cva_val = CVA_w_vars.S3_75thr.(f); % H x 2 
    cva_val_trans = cva_val'; % 
    S3_75thr_indicator.CVA_vars.weight.(f) = cva_val_trans; % 2 x H 
end
clear cva_val cva_val_trans
% mu, bands
% GMA, CVA lps
n_models = length(lp_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(lp_names{i}).S3_75thr.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end

for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);                         % n_models × H × B
    % GMA -----------------------------------------------------
    W = S3_75thr_indicator.GMA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_indicator.GMA_lps.irf.(f) = wirf_val;% H x B
    S3_75thr_indicator.GMA_lps.mu.(f) = mean(S3_75thr_indicator.GMA_lps.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.GMA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_indicator.GMA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
    % CVA ------------------------------------------------------
    W = S3_75thr_indicator.CVA_lps.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_indicator.CVA_lps.irf.(f) = wirf_val;% H x B
    S3_75thr_indicator.CVA_lps.mu.(f) = mean(S3_75thr_indicator.CVA_lps.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.CVA_lps.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_indicator.CVA_lps.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA vars
n_models = length(var_names);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names{i}).S3_75thr.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_75thr_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_indicator.GMA_vars.irf.(f) = wirf_val;% H x B
    S3_75thr_indicator.GMA_vars.mu.(f) = mean(S3_75thr_indicator.GMA_vars.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.GMA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_indicator.GMA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_75thr_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_indicator.CVA_vars.irf.(f) = wirf_val;
    S3_75thr_indicator.CVA_vars.mu.(f) = mean(S3_75thr_indicator.CVA_vars.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.CVA_vars.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_indicator.CVA_vars.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
% GMA, CVA girf
n_models = length(var_names_rob);
for j = 1:length(fields)
    f = fields{j};
    % without placeholder, it causes error in VAR version
    tmp = zeros(n_models,settings.est.IRF_hor,settings.est.n_btstrp);
    for i = 1:n_models
        tmp(i,:,:) = irfpaths.(var_names_rob{i}).S3_75thr.(f);
        % H x n_MC -> 1 x H x n_MC
    end
    irf_stack.(f) = tmp;
end
for j = 1:length(fields)
    f = fields{j};
    S = irf_stack.(f);
   % --------------------GMA---------------------------
    W = S3_75thr_indicator.GMA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_indicator.GMA_vars_girf.irf.(f) = wirf_val;% H x B
    S3_75thr_indicator.GMA_vars_girf.mu.(f) = mean(S3_75thr_indicator.GMA_vars_girf.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.GMA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_indicator.GMA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
   % --------------------CVA---------------------------
    W = S3_75thr_indicator.CVA_vars.weight.(f);    % n_models × H
    wirf_val = squeeze(sum(W.* S, 1));  % H x B
    S3_75thr_indicator.CVA_vars_girf.irf.(f) = wirf_val;
    S3_75thr_indicator.CVA_vars_girf.mu.(f) = mean(S3_75thr_indicator.CVA_vars_girf.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.CVA_vars_girf.band_lo.(f) = prctile(wirf_val,pct_lo,2);
    S3_75thr_indicator.CVA_vars_girf.band_hi.(f) = prctile(wirf_val,pct_hi,2);
end
%
clear wirf_val irf_stack n_models S W
% --------------------------------------------------------------------
% mavg: aLP 
% --------------------------------------------------------------------
% RS, MSPE alpha
% settings 
n_models_lp  = length(lp_names);
n_models_var = length(var_names);

rs_sumlp = struct;
mspe_sumlp =struct;
rs_sumvar = struct;
mspe_sumvar =struct;
%
for j = 1:length(fields) 
    f = fields{j};
    rs_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumlp.(f) = zeros(settings.est.IRF_hor, 1);
    % lps
    tss.S3_75thr.state0_vol_empl = tss.S3_75thr.state0_lps_empl;
    tss.S3_75thr.state0_ff_empl = tss.S3_75thr.state0_lps_empl;
    tss.S3_75thr.state1_vol_empl = tss.S3_75thr.state1_lps_empl;
    tss.S3_75thr.state1_ff_empl = tss.S3_75thr.state1_lps_empl;
    tss.S3_75thr.state0_vol_ip = tss.S3_75thr.state0_lps_ip;
    tss.S3_75thr.state0_ff_ip = tss.S3_75thr.state0_lps_ip;
    tss.S3_75thr.state1_vol_ip = tss.S3_75thr.state1_lps_ip;
    tss.S3_75thr.state1_ff_ip = tss.S3_75thr.state1_lps_ip;
    for i = 1:n_models_lp
        % RSQ
        rs_val = 1 - (ssr.(lp_names{i}).S3_75thr.(f) ./ tss.S3_75thr.(f)); % H x 1
        rs_sumlp.(f) = rs_sumlp.(f) + rs_val;
        % MSPE (horizon-specific)
        mspe_val = zeros(settings.est.IRF_hor, 1);
        for h = 1:settings.est.IRF_hor
            mspe_val(h,:) = (1/(T_62-settings.est.n_lag_short-h))*ssr.(lp_names{i}).S3_75thr.(f)(h,:); % H x 1
            mspe_sumlp.(f)(h,:) = mspe_sumlp.(f)(h,:) + mspe_val(h,:);
        end
    end
    % vars -------------------------------------------------------------
    rs_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    mspe_sumvar.(f) = zeros(settings.est.IRF_hor, 1);
    tss.S3_75thr.state0_vol_empl = tss.S3_75thr.state0_vars_empl;
    tss.S3_75thr.state0_ff_empl = tss.S3_75thr.state0_vars_empl;
    tss.S3_75thr.state1_vol_empl = tss.S3_75thr.state1_vars_empl;
    tss.S3_75thr.state1_ff_empl = tss.S3_75thr.state1_vars_empl;
    tss.S3_75thr.state0_vol_ip = tss.S3_75thr.state0_vars_ip;
    tss.S3_75thr.state0_ff_ip = tss.S3_75thr.state0_vars_ip;
    tss.S3_75thr.state1_vol_ip = tss.S3_75thr.state1_vars_ip;
    tss.S3_75thr.state1_ff_ip = tss.S3_75thr.state1_vars_ip;
    tss.S3_75thr.linear_vol_empl = tss.S3_75thr.linear_vars_empl;
    tss.S3_75thr.linear_vol_ip = tss.S3_75thr.linear_vars_ip;
    tss.S3_75thr.linear_ff_empl = tss.S3_75thr.linear_vars_empl;
    tss.S3_75thr.linear_ff_ip = tss.S3_75thr.linear_vars_ip;
    for i = 1:n_models_var
        rs_val = 1 - (ssr.(var_names{i}).S3_75thr.(f) ./ tss.S3_75thr.(f));
        rs_sumvar.(f) = rs_sumvar.(f) + rs_val;
        % MSPE (fixed denominator for VAR)
        mspe_val = (1/(T_62-settings.est.n_lag))*ssr.(var_names{i}).S3_75thr.(f); % H x 1
        mspe_sumvar.(f) = mspe_sumvar.(f) + mspe_val;
    end
end
% alpLP, weights
for j = 1:length(fields)
    f = fields{j};
    S3_75thr_indicator.aLP_RS.alpha.(f) = rs_sumlp.(f)./(rs_sumlp.(f) + rs_sumvar.(f)); % H x 1
    S3_75thr_indicator.aLP_MSPE.alpha.(f) = mspe_sumvar.(f)./(mspe_sumlp.(f) + mspe_sumvar.(f));
    % average weight for each horizon 
    % aLP_RS_GMA
    % Multiply: [H x 1] .* [2 x H ]' = [H x 2]
    weight_lps = S3_75thr_indicator.aLP_RS.alpha.(f) .* S3_75thr_indicator.GMA_lps.weight.(f)'; % H x 2
    alpha_var = 1 - S3_75thr_indicator.aLP_RS.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_75thr_indicator.GMA_vars.weight.(f)'; % H x 2
    S3_75thr_indicator.aLP_RS_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_RS_CVA
    weight_lps = S3_75thr_indicator.aLP_RS.alpha.(f) .* S3_75thr_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_75thr_indicator.CVA_vars.weight.(f)';
    S3_75thr_indicator.aLP_RS_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_GMA
    weight_lps = S3_75thr_indicator.aLP_MSPE.alpha.(f) .* S3_75thr_indicator.GMA_lps.weight.(f)';
    alpha_var = 1 - S3_75thr_indicator.aLP_MSPE.alpha.(f); % H x 1
    weight_vars = alpha_var .* S3_75thr_indicator.GMA_vars.weight.(f)';
    S3_75thr_indicator.aLP_MSPE_GMA.weight.(f) = [weight_vars weight_lps]; % H x 4
    % aLP_MSPE_CVA
    weight_lps = S3_75thr_indicator.aLP_MSPE.alpha.(f) .* S3_75thr_indicator.CVA_lps.weight.(f)';
    weight_vars = alpha_var .* S3_75thr_indicator.CVA_vars.weight.(f)';
    S3_75thr_indicator.aLP_MSPE_CVA.weight.(f) = [weight_vars weight_lps]; % H x 4
end
clear n_models_var n_models_lp weight_lps weight_vars alpha_var
clear mspe_sumvar mspe_sumlp mspe_val rs_sumvar rs_sumlp rs_val
% irf, mu, bands
for j = 1:length(fields)
    f = fields{j};
    S3_75thr_indicator.aLP_RS_GMA.irf.(f) = S3_75thr_indicator.aLP_RS.alpha.(f) .* S3_75thr_indicator.GMA_lps.irf.(f)...
    + (1-S3_75thr_indicator.aLP_RS.alpha.(f)) .* S3_75thr_indicator.GMA_vars.irf.(f); 
    S3_75thr_indicator.aLP_RS_GMA.mu.(f) = mean(S3_75thr_indicator.aLP_RS_GMA.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.aLP_RS_GMA.band_lo.(f) = prctile(S3_75thr_indicator.aLP_RS_GMA.irf.(f),pct_lo,2);
    S3_75thr_indicator.aLP_RS_GMA.band_hi.(f) = prctile(S3_75thr_indicator.aLP_RS_GMA.irf.(f),pct_hi,2);
    %
    S3_75thr_indicator.aLP_RS_CVA.irf.(f) = S3_75thr_indicator.aLP_RS.alpha.(f) .* S3_75thr_indicator.CVA_lps.irf.(f)...
    + (1-S3_75thr_indicator.aLP_RS.alpha.(f)) .* S3_75thr_indicator.CVA_vars.irf.(f); 
    S3_75thr_indicator.aLP_RS_CVA.mu.(f) = mean(S3_75thr_indicator.aLP_RS_CVA.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.aLP_RS_CVA.band_lo.(f) = prctile(S3_75thr_indicator.aLP_RS_CVA.irf.(f),pct_lo,2);
    S3_75thr_indicator.aLP_RS_CVA.band_hi.(f) = prctile(S3_75thr_indicator.aLP_RS_CVA.irf.(f),pct_hi,2);
    %
    S3_75thr_indicator.aLP_MSPE_GMA.irf.(f) = S3_75thr_indicator.aLP_MSPE.alpha.(f) .* S3_75thr_indicator.GMA_lps.irf.(f)...
    + (1-S3_75thr_indicator.aLP_MSPE.alpha.(f)) .* S3_75thr_indicator.GMA_vars.irf.(f); 
    S3_75thr_indicator.aLP_MSPE_GMA.mu.(f) = mean(S3_75thr_indicator.aLP_MSPE_GMA.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.aLP_MSPE_GMA.band_lo.(f) = prctile(S3_75thr_indicator.aLP_MSPE_GMA.irf.(f),pct_lo,2);
    S3_75thr_indicator.aLP_MSPE_GMA.band_hi.(f) = prctile(S3_75thr_indicator.aLP_MSPE_GMA.irf.(f),pct_hi,2);
    %
    S3_75thr_indicator.aLP_MSPE_CVA.irf.(f) = S3_75thr_indicator.aLP_MSPE.alpha.(f) .* S3_75thr_indicator.CVA_lps.irf.(f)...
    + (1-S3_75thr_indicator.aLP_MSPE.alpha.(f)) .* S3_75thr_indicator.CVA_vars.irf.(f); 
    S3_75thr_indicator.aLP_MSPE_CVA.mu.(f) = mean(S3_75thr_indicator.aLP_MSPE_CVA.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.aLP_MSPE_CVA.band_lo.(f) = prctile(S3_75thr_indicator.aLP_MSPE_CVA.irf.(f),pct_lo,2);
    S3_75thr_indicator.aLP_MSPE_CVA.band_hi.(f) = prctile(S3_75thr_indicator.aLP_MSPE_CVA.irf.(f),pct_hi,2);
    % --------------- Robustness: replace SVAR with GIRF-------------------
    S3_75thr_indicator.aLP_RS_GMA_girf.irf.(f) = S3_75thr_indicator.aLP_RS.alpha.(f) .* S3_75thr_indicator.GMA_lps.irf.(f)...
    + (1-S3_75thr_indicator.aLP_RS.alpha.(f)) .* S3_75thr_indicator.GMA_vars_girf.irf.(f); 
    S3_75thr_indicator.aLP_RS_GMA_girf.mu.(f) = mean(S3_75thr_indicator.aLP_RS_GMA_girf.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.aLP_RS_GMA_girf.band_lo.(f) = prctile(S3_75thr_indicator.aLP_RS_GMA_girf.irf.(f),pct_lo,2);
    S3_75thr_indicator.aLP_RS_GMA_girf.band_hi.(f) = prctile(S3_75thr_indicator.aLP_RS_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_75thr_indicator.aLP_RS_CVA_girf.irf.(f) = S3_75thr_indicator.aLP_RS.alpha.(f) .* S3_75thr_indicator.CVA_lps.irf.(f)...
    + (1-S3_75thr_indicator.aLP_RS.alpha.(f)) .* S3_75thr_indicator.CVA_vars_girf.irf.(f); 
    S3_75thr_indicator.aLP_RS_CVA_girf.mu.(f) = mean(S3_75thr_indicator.aLP_RS_CVA_girf.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.aLP_RS_CVA_girf.band_lo.(f) = prctile(S3_75thr_indicator.aLP_RS_CVA_girf.irf.(f),pct_lo,2);
    S3_75thr_indicator.aLP_RS_CVA_girf.band_hi.(f) = prctile(S3_75thr_indicator.aLP_RS_CVA_girf.irf.(f),pct_hi,2);
    %
    S3_75thr_indicator.aLP_MSPE_GMA_girf.irf.(f) = S3_75thr_indicator.aLP_MSPE.alpha.(f) .* S3_75thr_indicator.GMA_lps.irf.(f)...
    + (1-S3_75thr_indicator.aLP_MSPE.alpha.(f)) .* S3_75thr_indicator.GMA_vars_girf.irf.(f); 
    S3_75thr_indicator.aLP_MSPE_GMA_girf.mu.(f) = mean(S3_75thr_indicator.aLP_MSPE_GMA_girf.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.aLP_MSPE_GMA_girf.band_lo.(f) = prctile(S3_75thr_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_lo,2);
    S3_75thr_indicator.aLP_MSPE_GMA_girf.band_hi.(f) = prctile(S3_75thr_indicator.aLP_MSPE_GMA_girf.irf.(f),pct_hi,2);
    %
    S3_75thr_indicator.aLP_MSPE_CVA_girf.irf.(f) = S3_75thr_indicator.aLP_MSPE.alpha.(f) .* S3_75thr_indicator.CVA_lps.irf.(f)...
    + (1-S3_75thr_indicator.aLP_MSPE.alpha.(f)) .* S3_75thr_indicator.CVA_vars_girf.irf.(f); 
    S3_75thr_indicator.aLP_MSPE_CVA_girf.mu.(f) = mean(S3_75thr_indicator.aLP_MSPE_CVA_girf.irf.(f), 2, "omitmissing");
    S3_75thr_indicator.aLP_MSPE_CVA_girf.band_lo.(f) = prctile(S3_75thr_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_lo,2);
    S3_75thr_indicator.aLP_MSPE_CVA_girf.band_hi.(f) = prctile(S3_75thr_indicator.aLP_MSPE_CVA_girf.irf.(f),pct_hi,2);
end

%%
save("Combine_Results_Emp.mat",'S3_75thr_indicator','settings','-append')