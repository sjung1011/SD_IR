function [IRF_draws, IRF_hat] = ...
         GIRF_VAR_boot_stateDep(Y_raw, S_raw, settings)
% State-dependent GIRF.
% It follows weight of SVAR: they are substitutes (not used together) in this experiment
% ------------------------------------------------------------------
% 1. Settings & book-keeping
% ------------------------------------------------------------------
respV.empl = settings.est.resp_vars(1);
respV.ip  = settings.est.resp_vars(2);

recurShock.vol = settings.est.recurShock_vol;
recurShock.ff  = settings.est.recurShock_ff;
%normalizeV.vol = recurShock.vol;
%normalizeV.ff  = recurShock.ff;

%H  = settings.est.IRF_hor;
%B  = settings.est.n_btstrp;
%p  = settings.est.n_lag;

shocks = {'vol','ff'};
resps  = {'empl','ip'};

IRF_draws = struct();  IRF_hat = struct();

% ------------------------------------------------------------------
% 2. Loop over states, shocks, responses
% ------------------------------------------------------------------
for i_state = 0:1
    idx   = S_raw == i_state;
    Y_sub = Y_raw(idx,:);        % T × n  for this state
    
    for s = 1:numel(shocks)
        shockName  = shocks{s};
        shock_col  = recurShock.(shockName);
        %norm_col   = normalizeV.(shockName);

        for r = 1:numel(resps)
            respName = resps{r};
            resp_col = respV.(respName);

            field = ['state' num2str(i_state) '_' shockName '_' respName];

            % ---------- 2a. estimate ----------
            % delta is set following settings.est.delta_mode
            [irf_path,irf_hat, ~, ~] = GIRF_VAR_boot(Y_sub, settings,resp_col,shock_col);
            IRF_draws.(field) = irf_path;
            IRF_hat.(field) = irf_hat;
        end
    end
end
end