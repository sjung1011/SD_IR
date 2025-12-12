function [IRF, GIC, ei, SSR, SE] = BLP_est_stateDep(Y_raw, S_raw, settings)
% Compute state-dependent BLP IRFs for government spending (g) and output (y)
% in two regimes (S = 0, S = 1), normalized by the shock to x (1st col of Y).
% Inputs:
%   Y:       T x 3 data matrix [x g y]
%   S:       T x 1 binary state vector
%   settings: structure (should include .est.IRF_hor, .est.n_lag, .est.resp_vars, .est.shock)
% Outputs:
%   IRF_g, IRF_y: IRFs (IRF_hor x 2) for g and y, columns are [S==0, S==1]
%   extra_outputs: struct, can include diagnostics for each regime/var

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
        [IRF_g_tmp, GIC_g_tmp, ei_g_tmp, SSR_g_tmp, SE_g_tmp] = BLP_est(Y_sub, settings, respV.g, recurShock);
        IRF.(['state' num2str(i_state) '_g']) = IRF_g_tmp;
        GIC.(['state' num2str(i_state) '_g']) = GIC_g_tmp;
        ei.(['state' num2str(i_state) '_g']) = ei_g_tmp;
        SSR.(['state' num2str(i_state) '_g']) = SSR_g_tmp;
        SE.(['state' num2str(i_state) '_g']) = SE_g_tmp';
        % y
        [IRF_y_tmp, GIC_y_tmp, ei_y_tmp, SSR_y_tmp, SE_y_tmp] = BLP_est(Y_sub, settings, respV.y, recurShock);
        IRF.(['state' num2str(i_state) '_y']) = IRF_y_tmp;
        GIC.(['state' num2str(i_state) '_y']) = GIC_y_tmp;
        ei.(['state' num2str(i_state) '_y']) = ei_y_tmp;
        SSR.(['state' num2str(i_state) '_y']) = SSR_y_tmp;
        SE.(['state' num2str(i_state) '_y']) = SE_y_tmp';
    end
end
