function [IRF, BIC, ei, SSR, TSS, SE] = BVAR_est_stateDep(Y_raw, S_raw, settings)
    % settings 
    respV.g = settings.est.resp_vars(1); %2 in Y_raw
    respV.y = settings.est.resp_vars(2); %3
    recurShock = settings.est.shock; 
    normalizeV = recurShock;
    IRF_hor = settings.est.IRF_hor;
    
    % bias_corrected =1 if BC
    % Preallocate
    IRF =struct;
    BIC =struct;
    ei =struct;
    SSR =struct;
    TSS =struct;
    SE =struct;
    

    for i_state = 0:1
        idx = (S_raw == i_state);
        Y_sub = Y_raw(idx,:);
        % g
        [IRF_g_tmp, BIC_g_tmp, ei_g_tmp, SSR_g_tmp,SE_g_tmp,~] ...
            = BVAR_est(Y_sub, settings, respV.g,recurShock,normalizeV);
        IRF.(['state' num2str(i_state) '_g']) = IRF_g_tmp;
        BIC.(['state' num2str(i_state) '_g']) = BIC_g_tmp;
        ei.(['state' num2str(i_state) '_g']) = ei_g_tmp;
        SSR.(['state' num2str(i_state) '_g']) = SSR_g_tmp;
        SE.(['state' num2str(i_state) '_g']) = SE_g_tmp';
        % y
        [IRF_y_tmp, BIC_y_tmp, ei_y_tmp, SSR_y_tmp,SE_y_tmp,~] ...
            = BVAR_est(Y_sub, settings, respV.y,recurShock,normalizeV);
        IRF.(['state' num2str(i_state) '_y']) = IRF_y_tmp;
        BIC.(['state' num2str(i_state) '_y']) = BIC_y_tmp;
        ei.(['state' num2str(i_state) '_y']) = ei_y_tmp;
        SSR.(['state' num2str(i_state) '_y']) = SSR_y_tmp;
        SE.(['state' num2str(i_state) '_y']) = SE_y_tmp';
        
    end
end