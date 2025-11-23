function [CVA_w] = CVA_est_lp_linear_emp(Y_raw, settings,ei_lp,ei_lp_contmp)
    % settings 
    %respV.empl = settings.est.resp_vars(1);
    %respV.ip  = settings.est.resp_vars(2);
    recurShock.vol = settings.est.recurShock_vol;
    recurShock.ff  = settings.est.recurShock_ff;
    
    % Preallocate
    shocks = {'vol','ff'};
    resps  = {'empl','ip'};
    CVA_w = struct();
 % ------------------------------------------------------------------
 % 2. Loop over states, shocks, responses
 % ------------------------------------------------------------------

    
    for s = 1:numel(shocks)
        shockName  = shocks{s};
        shock_col  = recurShock.(shockName);
        
        for r = 1:numel(resps)
            respName = resps{r};
            %resp_col = respV.(respName);
            field = ['linear_'  shockName '_' respName];
            % ---------- estimate ----------
            [CVA_w_tmp] = CVA_est_lp_emp(Y_raw,settings,ei_lp.(field),ei_lp_contmp.(field),shock_col);
            CVA_w.(field) = CVA_w_tmp';
        end
    end
end

