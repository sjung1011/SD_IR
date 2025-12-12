function [IRF, BIC_bvar, Res, SSR, SE_B, IRF_draws] = BVAR_est(Y_raw,settings, responseV,recurShock,normalizeV);
% Function for estimating IRFs using a Bayesian VAR approach

% preparations
IRF_hor = settings.est.IRF_hor;
nlags = settings.est.n_lag;
prior = settings.est.prior;
ndraw = settings.est.posterior_ndraw;


% estimate Bayesian VAR and get posterior mean and variance for VAR coef.

if ndraw == 0 % only use posterior mean of VAR coef
    [~,By,Sigma] = BVAR(Y_raw,nlags,prior);
else % use posterior mean and variance of VAR coef
    [~,~,Sigma,~,VAR_coef_post_mean,VAR_coef_post_vce_inv] = BVAR(Y_raw,nlags,prior);
end
G = chol(Sigma, 'lower'); % Warning: correspond to matrix C in our paper
ShockVector = G(:,recurShock);

% estimate IRF

if ndraw == 0 % only use posterior mean of VAR coef to compute IRF
    IRF = IRF_SVAR(By,ShockVector,IRF_hor - 1); % IRF to one unit of shock
    IRF_draws = zeres(IRF_hor,settings.est.posterior_ndraw);
else % compute posterior mean of IRF based on posterior draws
    IRF_draws = IRF_BVAR(VAR_coef_post_mean,VAR_coef_post_vce_inv,ShockVector,IRF_hor - 1,ndraw,Y_raw); % IRF draws
    IRF = mean(IRF_draws,3); % posterior mean of IRF draws
    IRF_draws = IRF_draws(responseV,:,:);
    IRF_draws = squeeze(IRF_draws);
end

IRF = IRF(responseV,:) / IRF(normalizeV,1); % normalize by response of normalization variable
IRF = IRF';
% BIC if ndraws>0 (we didn't set code of ndraw==0)
[~,~,~,~,VAR_coef_post_mean,VAR_coef_post_vce_inv] = BVAR(Y_raw,nlags,prior);
[~, BIC_bvar,Res, SE_B] = IRF_BVAR(VAR_coef_post_mean,VAR_coef_post_vce_inv,ShockVector,IRF_hor - 1,ndraw,Y_raw);
% for alpha LP guidance 
Res_sq = Res.^2; % Square each element in Res
SSR_all = sum(Res_sq,1)'; % (1 x n)' sum residuals across t 
SSR_mat = repmat(SSR_all,1,IRF_hor); % n x h, copy sumRes as VARs has horizon-specific Res
SSR = SSR_mat(responseV,:); % 1 x h
SSR = SSR'; % h x 1, stored by model
% error band in empirical: asymptotic variance of VAR coeff
SE_B = SE_B(1+recurShock,responseV); % 1 x 1
SE_B = repmat(SE_B,1,IRF_hor); % 1 x h
end