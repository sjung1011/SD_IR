function [IRF, band_lo, band_hi, ei] = GIRF_VAR_est_stateDep(Y_raw, S_raw, settings)
    % settings 
    respV.g = settings.est.resp_vars(1); %2 in Y_raw
    respV.y = settings.est.resp_vars(2); %3
    recurShock = settings.est.shock; 
    
    % Preallocate
    IRF =struct;
    band_lo =struct;
    band_hi =struct;
    ei =struct;

    for i_state = 0:1
        idx = (S_raw == i_state);
        Y_sub = Y_raw(idx,:);
        % g
        [IRF_g_tmp, band_lo_g_tmp, band_hi_g_tmp, ei_g_tmp] ...
            = GIRF_VAR_est(Y_sub, settings, respV.g,recurShock,1);
        IRF.(['state' num2str(i_state) '_g']) = IRF_g_tmp;
        band_lo.(['state' num2str(i_state) '_g']) = band_lo_g_tmp;
        band_hi.(['state' num2str(i_state) '_g']) = band_hi_g_tmp;
        ei.(['state' num2str(i_state) '_g']) = ei_g_tmp;
        % y
        [IRF_y_tmp, band_lo_y_tmp, band_hi_y_tmp, ei_y_tmp] ...
            = GIRF_VAR_est(Y_sub, settings, respV.y,recurShock,1);
        IRF.(['state' num2str(i_state) '_y']) = IRF_y_tmp;
        band_lo.(['state' num2str(i_state) '_y']) = band_lo_y_tmp;
        band_hi.(['state' num2str(i_state) '_y']) = band_hi_y_tmp;
        ei.(['state' num2str(i_state) '_y']) = ei_y_tmp;
    end
end