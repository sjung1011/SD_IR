function [IRF, GIC, ei, SSR, SE] ...
    = LP_shrink_est_stateDep(Y_raw, S_raw, settings)
% Computes state-dependent IRFs for g and y using penalized LP for each regime.
% Y_raw: T x 3 (columns: x, g, y)
% S_raw: T x 1 (state dummy, 0/1)
% settings: struct with estimation settings
    % settings 
    respV.g = settings.est.resp_vars(1); %2 in Y_raw
    respV.y = settings.est.resp_vars(2); %3
    recurShock = settings.est.shock; 
    normalizeV = recurShock;
    IRF_hor = settings.est.IRF_hor;

    % Preallocate
    IRF =struct;
    GIC =struct;
    ei =struct;
    SSR =struct;
    SE =struct;
    

    for i_state = 0:1
        idx = (S_raw == i_state);
        Y_sub = Y_raw(idx,:);
        % g
        [IRF_g_tmp, GIC_g_tmp, ei_g_tmp, SSR_g_tmp, SE_g_tmp]...
            = LP_shrink_est(Y_sub, settings, respV.g,recurShock, normalizeV);
        IRF.(['state' num2str(i_state) '_g']) = IRF_g_tmp;
        GIC.(['state' num2str(i_state) '_g']) = GIC_g_tmp;
        ei.(['state' num2str(i_state) '_g']) = ei_g_tmp;
        SSR.(['state' num2str(i_state) '_g']) = SSR_g_tmp;
        SE.(['state' num2str(i_state) '_g']) = SE_g_tmp';
        % y
        [IRF_y_tmp, GIC_y_tmp, ei_y_tmp, SSR_y_tmp, SE_y_tmp] ...
            = LP_shrink_est(Y_sub, settings, respV.y,recurShock, normalizeV);
        IRF.(['state' num2str(i_state) '_y']) = IRF_y_tmp;
        GIC.(['state' num2str(i_state) '_y']) = GIC_y_tmp;
        ei.(['state' num2str(i_state) '_y']) = ei_y_tmp;
        SSR.(['state' num2str(i_state) '_y']) = SSR_y_tmp;
        SE.(['state' num2str(i_state) '_y']) = SE_y_tmp';
        
    end
    
end
