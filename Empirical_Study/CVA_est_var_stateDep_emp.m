function [CVA_w] = CVA_est_var_stateDep_emp(Y_raw, S_raw, settings, ei_var,ei_bvar)
    % settings 
    respV.empl = settings.est.resp_vars(1);
    respV.ip  = settings.est.resp_vars(2);
    
    % Preallocate
    shocks = {'vol','ff'};
    resps  = {'empl','ip'};
    CVA_w = struct();
 % ------------------------------------------------------------------
 % 2. Loop over states, shocks, responses
 % ------------------------------------------------------------------
for i_state = 0:1
    idx   = S_raw == i_state;
    Y_sub = Y_raw(idx,:);        % T × n  for this state
    
    for s = 1:numel(shocks)
        shockName  = shocks{s};
        
        for r = 1:numel(resps)
            respName = resps{r};
            resp_col = respV.(respName);
            field = ['state' num2str(i_state) '_' shockName '_' respName];

            % ---------- estimate ----------
            [CVA_w_tmp] = CVA_est_var_emp(Y_sub,settings,ei_var.(field),ei_bvar.(field),resp_col);
            CVA_w.(field) = CVA_w_tmp';
        end
    end
end
end

