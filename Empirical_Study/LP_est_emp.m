function [IRF_draws,IRF,GIC,Res,SSR,TSS] ...
= LP_est_emp(Y,settings,responseV,recurShock,normalizeV)
% Function for estimating IRFs using least-squares LP

% preparations
IRF_hor = settings.est.IRF_hor;
nlags = settings.est.n_lag_short; % short for LPs

% estimate IRF via LP
[IRF_resp,IRF_draws,nT_GIC,nv_GIC,est_Var_HAC,Res,SSR,TSS]... 
= IRF_LP_emp(Y,recurShock,responseV,nlags,IRF_hor - 1); % IRF to one unit of shock

% normalize
IRF_normalize = IRF_LP(Y,recurShock,normalizeV,nlags,0);
IRF = IRF_resp / IRF_normalize; % normalize by response of normalization variable
IRF = IRF';

%--------------------------------------
% GMA: Only compute GIC_ls here 
%--------------------------------------
% nT_BIC, nv_BIC, Res vary in h in each lag folder (thus, not reflect varying lags)
GIC = zeros(1,IRF_hor);

% go thru horizons
for h=1:IRF_hor
    GIC(:, h) = log(det(est_Var_HAC(h))) + nv_GIC(h) * log(nT_GIC(h)) / nT_GIC(h); 
end

