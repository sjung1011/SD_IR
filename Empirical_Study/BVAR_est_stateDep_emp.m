function [IRF_draws, IRF_hat, BIC, Res, SSR, SE] = ...
         BVAR_est_stateDep_emp(Y_raw, S_raw, settings)
% State-dependent BVAR with increased posterior draws.

% ------------------------------------------------------------------
% 1. Settings & book-keeping
% ------------------------------------------------------------------
respV.empl = settings.est.resp_vars(1);
respV.ip  = settings.est.resp_vars(2);

recurShock.vol = settings.est.recurShock_vol;
recurShock.ff  = settings.est.recurShock_ff;
normalizeV.vol = recurShock.vol;
normalizeV.ff  = recurShock.ff;

%H  = settings.est.IRF_hor;
%B  = settings.est.n_btstrp;
%p  = settings.est.n_lag;

shocks = {'vol','ff'};
resps  = {'empl','ip'};

IRF_draws = struct();  IRF_hat = struct();
BIC = struct();  Res = struct();
SSR = struct();  SE = struct();

% ------------------------------------------------------------------
% 2. Loop over states, shocks, responses
% ------------------------------------------------------------------
for i_state = 0:1
    idx   = S_raw == i_state;
    Y_sub = Y_raw(idx,:);        % T × n  for this state
    
    for s = 1:numel(shocks)
        shockName  = shocks{s};
        shock_col  = recurShock.(shockName);
        norm_col   = normalizeV.(shockName);

        for r = 1:numel(resps)
            respName = resps{r};
            resp_col = respV.(respName);

            field = ['state' num2str(i_state) '_' shockName '_' respName];

            % ---------- 2a. point estimate ----------
            [irf_hat, bic, Res_tmp, ssr, se, irf_draws] = ...
                BVAR_est(Y_sub, settings, resp_col, shock_col, norm_col);

            IRF_hat.(field) = irf_hat;
            BIC.(field)     = bic';
            Res.(field)     = Res_tmp;          % keep residuals if you wish
            SSR.(field)     = ssr;
            SE.(field)      = se';
            IRF_draws.(field) = irf_draws;
        end
    end
end
end