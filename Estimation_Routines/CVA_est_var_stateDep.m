function [CVA_w] = CVA_est_var_stateDep(Y_raw, S_raw, settings, ei_var,ei_bcvar,ei_bvar)
    % settings 
    respV.g = settings.est.resp_vars(1); %2 in Y_raw
    respV.y = settings.est.resp_vars(2); %3
    %recurShock = settings.est.shock; 
    %normalizeV = recurShock;
    %IRF_hor = settings.est.IRF_hor;
    
    
    % Preallocate
    CVA_w =struct;
    
        % g,y
        % state0
        idx = (S_raw == 0);
        Y_sub = Y_raw(idx,:);
        [CVA_w_g_tmp] = CVA_est_var(Y_sub,settings,ei_var.state0_g,ei_bcvar.state0_g,ei_bvar.state0_g,respV.g);
        CVA_w.state0_g = CVA_w_g_tmp'; % 4 x H
        [CVA_w_y_tmp] = CVA_est_var(Y_sub,settings,ei_var.state0_y,ei_bcvar.state0_y,ei_bvar.state0_y,respV.y);
        CVA_w.state0_y = CVA_w_y_tmp';
        % state1
        idx = (S_raw == 1);
        Y_sub = Y_raw(idx,:);
        [CVA_w_g_tmp] = CVA_est_var(Y_sub,settings,ei_var.state1_g,ei_bcvar.state1_g,ei_bvar.state1_g,respV.g);
        CVA_w.state1_g = CVA_w_g_tmp';
        [CVA_w_y_tmp] = CVA_est_var(Y_sub,settings,ei_var.state1_y,ei_bcvar.state1_y,ei_bvar.state1_y,respV.y);
        CVA_w.state1_y = CVA_w_y_tmp';
   
end