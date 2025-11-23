function [IRF_draws, IRF_hat, GIC, Res, SSR, TSS] = ...
         LP_est_linear_emp(Y_raw, settings)
% State-dependent BVAR with increased posterior draws.

% ------------------------------------------------------------------
% 1. Settings & book-keeping
% ------------------------------------------------------------------
respV.empl = settings.est.resp_vars(1);
respV.ip  = settings.est.resp_vars(2);

recurShock.vol = settings.est.recurShock_vol;
recurShock.ff  = settings.est.recurShock_ff;
normalizeV.vol = recurShock.vol;
normalizeV.ff  = settings.est.normalizeV_ff;

shocks = {'vol','ff'};
resps  = {'empl','ip'};

IRF_draws = struct();  IRF_hat = struct();
GIC = struct();  Res = struct();
SSR = struct();  TSS = struct();

% ------------------------------------------------------------------
% 2. Loop over states, shocks, responses
% ------------------------------------------------------------------
    
    for s = 1:numel(shocks)
        shockName  = shocks{s};
        shock_col  = recurShock.(shockName);
        norm_col   = normalizeV.(shockName);

        for r = 1:numel(resps)
            respName = resps{r};
            resp_col = respV.(respName);

            field = ['linear_' shockName '_' respName];

            % ---------- point estimate ----------
            [irf_draws,irf_hat,gic,res,ssr,tss]... 
            = LP_est_emp(Y_raw, settings, resp_col, shock_col, norm_col);
            
            IRF_hat.(field) = irf_hat;
            GIC.(field)     = gic';
            Res.(field)     = res;          % keep residuals if you wish
            SSR.(field)     = ssr;
            IRF_draws.(field) = irf_draws;
            TSS.(field)     = tss;
        end
    end
end
