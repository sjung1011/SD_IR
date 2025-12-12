function [IRF, GIC, ei, SSR, SE_HAC, lambda_opt] = ...
    LP_shrink_est(Y,settings ,responseV,recurShock,normalizeV)

% preparations
lambdaRange = settings.est.lambdaRange;
irfLimitOrder = settings.est.irfLimitOrder;
IRF_hor = settings.est.IRF_hor;
nlags = settings.est.n_lag;

% data for LP routine
[y, x, w] = LP_gen_data(Y,recurShock,responseV,nlags,0); % Warning: here w include both contemperaneous and lagged controls

% settings for LP shrinkage method
H_min = 0; % min horizon in IRF
H_max = IRF_hor - 1; % max horizon in IRF
r = irfLimitOrder + 1; % order of finite difference operator in the penalty term

% leave-one-out cross validation
nT = size(Y,1);
lambdaRange = lambdaRange * nT; % scale the grid of lambda (penalty strength) by nT
lambdaRange = [1e-4, lambdaRange, 1e10]; % allow regular OLS or completely smoothed
rss_cv = locproj_cv(y, x, w, H_min, H_max, r, lambdaRange, settings.est.CV_folds); % cross-validated MSE for each value of lambda
[~,lambda_opt_loc] = min(rss_cv);
lambda_opt = lambdaRange(lambda_opt_loc); % optimally tuned lambda

% re-estimate IRF via penalized LP using the full sample
[IRF_resp, ~, ~, ~, ~, nv_GIC,ei] = locproj(y, x, w, H_min, H_max, r, lambda_opt); % IRF to one unit of shock

IRF_normalize = IRF_LP(Y,recurShock,normalizeV,nlags,0); % response of normalization variable estimated by least-squares LP
IRF = IRF_resp / IRF_normalize; % normalize by response of normalization variable


%-----------------------------
% for GIC
%-----------------------------
GIC = zeros(1,IRF_hor); 
SSR = zeros(IRF_hor,1);
nT_GIC = zeros(1,IRF_hor);
nv_GIC = ones(1,IRF_hor)*nv_GIC;
est_Var_HAC =zeros(1,IRF_hor);
SE_HAC = zeros(1,IRF_hor);


% go thru horizons
for h = 0:IRF_hor-1
    field_name = sprintf('h%d', h); % Dynamically create the field name
    Res_sq=ei.(field_name).^2;
    SSR(h + 1, :) = sum(Res_sq, 1);
    nT_GIC(h+1) = size(ei.(field_name), 1); % Access the field and get the size
    [est_Var_HAC(:,h+1), SE_HAC(:,h+1)] = hac_standard_error(ei.(field_name),1,nT_GIC(h+1));
    GIC(:, h+1) = log(det(est_Var_HAC(:,h+1))) + nv_GIC(h+1) * log(nT_GIC(h+1)) / nT_GIC(h+1); % det is a scalar 
end
end