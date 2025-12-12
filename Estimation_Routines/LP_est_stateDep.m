function [IRF, GIC, ei, SSR, TSS, SE] = LP_est_stateDep(Y_raw, S_raw, settings, bias_corrected)
    % settings 
    respV.g = settings.est.resp_vars(1); %2 in Y_raw
    respV.y = settings.est.resp_vars(2); %3
    recurShock = settings.est.shock; 
    normalizeV = recurShock;
    IRF_hor = settings.est.IRF_hor;
    
    % bias_corrected =1 if BC
    % Preallocate
    %IRF_g = NaN(IRF_hor,2);
    %IRF_y = NaN(IRF_hor,2);
    % outputs = struct; % H x 1 x M after MCs
    IRF =struct;
    GIC =struct;
    ei =struct;
    SSR =struct;
    TSS =struct;
    SE =struct;
    

    for i_state = 0:1
        idx = (S_raw == i_state);
        Y_sub = Y_raw(idx,:);
        % g
        [IRF_g_tmp, GIC_g_tmp, ei_g_tmp, SSR_g_tmp,TSS_g_tmp, SE_g_tmp]...
            = LP_est(Y_sub, settings, respV.g,recurShock, normalizeV, bias_corrected);
        %IRF_g(:,i_state+1) = IRF_g_tmp;
        IRF.(['state' num2str(i_state) '_g']) = IRF_g_tmp;
        GIC.(['state' num2str(i_state) '_g']) = GIC_g_tmp;
        ei.(['state' num2str(i_state) '_g']) = ei_g_tmp;
        SSR.(['state' num2str(i_state) '_g']) = SSR_g_tmp;
        TSS.(['state' num2str(i_state) '_g']) = TSS_g_tmp;
        SE.(['state' num2str(i_state) '_g']) = SE_g_tmp';
        % y
        [IRF_y_tmp, GIC_y_tmp, ei_y_tmp, SSR_y_tmp,TSS_y_tmp, SE_y_tmp] ...
            = LP_est(Y_sub, settings, respV.y,recurShock, normalizeV, bias_corrected);
        %IRF_y(:,i_state+1) = IRF_y_tmp;
        IRF.(['state' num2str(i_state) '_y']) = IRF_y_tmp;
        GIC.(['state' num2str(i_state) '_y']) = GIC_y_tmp;
        ei.(['state' num2str(i_state) '_y']) = ei_y_tmp;
        SSR.(['state' num2str(i_state) '_y']) = SSR_y_tmp;
        TSS.(['state' num2str(i_state) '_y']) = TSS_y_tmp;
        SE.(['state' num2str(i_state) '_y']) = SE_y_tmp';
        
    end
end