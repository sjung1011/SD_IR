function [IRF,GIC,Res,SSR,TSS,SE_HAC] = LP_est(Y,settings,responseV,recurShock,normalizeV,bias_corrected)
% Function for estimating IRFs using least-squares LP

% preparations
IRF_hor = settings.est.IRF_hor;
nlags = settings.est.n_lag;

% estimate IRF via LP

[IRF_resp, w, nT_GIC, nv_GIC, est_Var_HAC, est_Var_HAC_bc...
    ,ei_ls, ei_bc,SSR_ls,SSR_bc, TSS, SE, SE_bc] = IRF_LP(Y,recurShock,responseV,nlags,IRF_hor - 1); % IRF to one unit of shock
Res = ei_ls;
SSR = SSR_ls;
SE_HAC = SE;
if bias_corrected == 1
    Res = ei_bc;
    SSR = SSR_bc;
    SE_HAC = SE_bc;
    [IRF_resp,~] = LP_CorrectBias(IRF_resp, w); % Herbst-Johanssen bias correction
end

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
    if bias_corrected == 1
        GIC(:, h) = log(det(est_Var_HAC_bc(h))) + nv_GIC(h) * log(nT_GIC(h)) / nT_GIC(h); 
    end
end

