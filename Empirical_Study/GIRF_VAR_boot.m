function [irf_boot, irf_mean, band_lo, band_hi] = GIRF_VAR_boot(Y, settings, responseV, recurShock, delta)
% GIRF_VAR_boot: Bootstrap confidence bands for GIRFs
%
% Inputs:
%   Y          (T × n) data
%   settings   struct with fields:
%                 .est.n_lag      = VAR order
%                 .est.IRF_hor    = horizon length (H+1)
%                 .est.n_boot     = # bootstrap replications
%   responseV  index of response variable
%   recurShock index of shocked variable
%   delta      shock size
%
% Outputs:
%   irf_boot   (H+1 × nBoot) bootstrap IRFs
%   irf_mean   (H+1 × 1)     mean across bootstrap IRFs
%   band_lo    (H+1 × 1)     5th percentile
%   band_hi    (H+1 × 1)     95th percentile
%

    nlags = settings.est.n_lag;
    H     = settings.est.IRF_hor - 1;
    nBoot = settings.est.n_btstrp;

    irf_boot = zeros(H+1, nBoot);

    % ---------- Estimate VAR once to get residuals ----------
    [Bc,By,~,~,Res,~,~,~,~] = VAR(Y, nlags);
    [T,n] = size(Y);
%
if nargin < 5 || isempty(delta)
    switch lower(settings.est.delta_mode)         % choose in settings
        case 'unit'          % structural IRF, Δ = 1
            delta = 1;
        case 'onesd'         % 1 s.d. of reduced-form residual
            delta = std(Res(:,recurShock));
        case 'twosd'         % 2 s.d. (non-linearities check)
            delta = 2*std(Res(:,recurShock));
    end
else
    delta = delta_in;        % user supplied in the call
end

% Wild bootstrap
for b = 1:nBoot
        % Wild multipliers: m_t ∈ {−1, +1} with prob 1/2
        m = 2*(rand(T-nlags,1) > 0.5) - 1;       % (T-nlags) × 1
        Res_b = Res .* m;                        % implicit expansion to (T-nlags) × n

        % Generate bootstrap sample Y_b using original coefficients and Res_b
        Y_b = bootstrap_VAR_sim(Bc, By, Y, Res_b, nlags);

        % Re-estimate VAR on bootstrap sample to capture parameter uncertainty
        [Bc_b,By_b,~,~,Res_hat_b,~,~,~,~] = VAR(Y_b, nlags);

        % Compute GIRF deterministically from re-estimated VAR and residuals
        irf_b = GIRF_from_VAR(By_b, Bc_b, Res_hat_b, responseV, recurShock, H, delta);
        irf_boot(:,b) = irf_b;
end
%{
% iid bootstrap
    for b = 1:nBoot
        % Step 1: bootstrap residuals (you can change to wild)
        idx = randi(size(Res,1), T-nlags, 1);
        Res_b = Res(idx,:);

        % Step 2: generate bootstrap sample Y_b
        Y_b = bootstrap_VAR_sim(Bc,By,Y,Res_b,nlags);

        % Step 3: re-estimate VAR on bootstrap sample
        [Bc_b,By_b,~,~,~,~,~,~,~] = VAR(Y_b,nlags);

        % Step 4: compute GIRF for this bootstrap draw
        irf_b = GIRF_from_VAR(By_b,Bc_b,Res_b,responseV,recurShock,H,delta);
        irf_boot(:,b) = irf_b;
    end
%}
    % ---------- Summarize ----------
    irf_mean = mean(irf_boot,2);
    band_lo  = prctile(irf_boot,5,2);
    band_hi  = prctile(irf_boot,95,2);
end
