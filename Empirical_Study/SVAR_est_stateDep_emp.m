function [IRF_boot, IRF_hat, BIC, Res, SSR, TSS] = ...
         SVAR_est_stateDep_emp(Y_raw, S_raw, settings, bias_corrected)
% State-dependent SVAR with wild (Rademacher) bootstrap IRFs.
% Uses the residuals produced by SVAR_est.

% ------------------------------------------------------------------
% 1. Settings & book-keeping
% ------------------------------------------------------------------
respV.empl = settings.est.resp_vars(1);
respV.ip  = settings.est.resp_vars(2);

recurShock.vol = settings.est.recurShock_vol;
recurShock.ff  = settings.est.recurShock_ff;
normalizeV.vol = settings.est.normalizeV_vol;
normalizeV.ff  = settings.est.normalizeV_ff;

H  = settings.est.IRF_hor;
B  = settings.est.n_btstrp;
p  = settings.est.n_lag;

shocks = {'vol','ff'};
resps  = {'empl','ip'};

IRF_boot = struct();  IRF_hat = struct();
BIC = struct();  Res = struct();
SSR = struct();  TSS = struct();  %SE = struct();

% ------------------------------------------------------------------
% 2. Loop over states, shocks, responses
% ------------------------------------------------------------------
for i_state = 0:1
    idx   = S_raw == i_state;
    Y_sub = Y_raw(idx,:);        % T × n  for this state
    T     = size(Y_sub,1);

    for s = 1:numel(shocks)
        shockName  = shocks{s};
        shock_col  = recurShock.(shockName);
        norm_col   = normalizeV.(shockName);

        for r = 1:numel(resps)
            respName = resps{r};
            resp_col = respV.(respName);

            field = ['state' num2str(i_state) '_' shockName '_' respName];

            % ---------- 2a. point estimate ----------
            [irf_hat, bic, Res_tmp, ssr, tss, ~] = ...
                SVAR_est(Y_sub, settings, resp_col, shock_col, norm_col, bias_corrected);

            IRF_hat.(field) = irf_hat;
            BIC.(field)     = bic';
            Res.(field)     = Res_tmp;          % keep residuals if you wish
            SSR.(field)     = ssr;
            TSS.(field)     = tss;
            %SE.(field)      = se';

            % ---------- 2b. wild bootstrap ----------
            IRF_boot.(field) = NaN(H, B);

            for b = 1:B
                Yb = wild_bootstrap_VAR(Y_sub, Res_tmp, p, T);
                irf_b = SVAR_est(Yb, settings, resp_col, shock_col, norm_col, bias_corrected);
                IRF_boot.(field)(:,b) = irf_b;
            end
        end
    end
end
end