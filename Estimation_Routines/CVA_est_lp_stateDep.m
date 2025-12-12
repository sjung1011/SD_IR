function [CVA_w] = CVA_est_lp_stateDep(Y_raw, S_raw, settings, ei_lp,ei_bclp,ei_blp,ei_pen)
    % settings 
    respV.g = settings.est.resp_vars(1); %2 in Y_raw
    respV.y = settings.est.resp_vars(2); %3
    recurShock = settings.est.shock; 
    %normalizeV = recurShock;
    %IRF_hor = settings.est.IRF_hor;
    
    
    % Preallocate
    CVA_w =struct;
    
        % g,y
        % state0
        idx = (S_raw == 0);
        Y_sub = Y_raw(idx,:);
        [CVA_w_g_tmp] = CVA_est_lp(Y_sub,settings,ei_lp.state0_g,ei_bclp.state0_g,ei_blp.state0_g,ei_pen.state0_g,respV.g,recurShock);
        CVA_w.state0_g = CVA_w_g_tmp'; % 4 x H
        [CVA_w_y_tmp] = CVA_est_lp(Y_sub,settings,ei_lp.state0_y,ei_bclp.state0_y,ei_blp.state0_y,ei_pen.state0_y,respV.y,recurShock);
        CVA_w.state0_y = CVA_w_y_tmp';
        % state1
        idx = (S_raw == 1);
        Y_sub = Y_raw(idx,:);
        [CVA_w_g_tmp] = CVA_est_lp(Y_sub,settings,ei_lp.state1_g,ei_bclp.state1_g,ei_blp.state1_g,ei_pen.state1_g,respV.g,recurShock);
        CVA_w.state1_g = CVA_w_g_tmp';
        [CVA_w_y_tmp] = CVA_est_lp(Y_sub,settings,ei_lp.state1_y,ei_bclp.state1_y,ei_blp.state1_y,ei_pen.state1_y,respV.y,recurShock);
        CVA_w.state1_y = CVA_w_y_tmp';
   
end