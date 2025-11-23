function params = estimateParams_RZ(filename, sample_size, dgp_type, data_trans)
    % Estimate parameters for state-dependent SVAR(4) model from real data
    %
    % Inputs:
    % filename     : path to RZDAT.xlsx
    % sample_size  : "long"  (1889q1–2015q4)  or "short" (1947q1–2015q4)
    % dgp_type     : "TVAR", "MSVAR", or "STVAR"
    % Outputs:
    % params - SVAR (C) and reduced-form VAR parameters (A,k)

    % ---------- Load data and construct series
    sheetname = 'rzdat';
    TBL = readtable(filename,'Sheet',sheetname,'PreserveVariableNames',true);
    % ----------- Choose data transformation
    if data_trans == "level"
        pgdpLag  = lagmatrix(TBL.pgdp,1);
        ynormLag = lagmatrix(TBL.rgdp_pott6,1);
        x = TBL.news ./ (ynormLag .* pgdpLag);  % military news
        g = log(TBL.ngov ./ TBL.pgdp);   % Log of real government spending in level
         y = log(TBL.rgdp);   %  in level: doens't work due to trend
    else % stationary
        pgdpLag  = lagmatrix(TBL.pgdp,1);
        ynormLag = lagmatrix(TBL.rgdp_pott6,1);
        x = TBL.news ./ (ynormLag .* pgdpLag);                     % military news
        g = TBL.ngov ./ (TBL.pgdp .* TBL.rgdp_pott6);              % government spending
        y = TBL.rgdp ./ TBL.rgdp_pott6;                       % real GDP
    end
    % ----------Choose sample
    sample_size = string(sample_size);    

    if sample_size == "long" % 1889q1 (1890q1, drop NaN)–2015q4
       idx = find(~isnan(x));               % logical index over the whole vector
    else % 1947q1–2015q4 (289:564)
        rawIdx = 289:find(~isnan(x),1,'last');
        idx    = rawIdx(~isnan(x(rawIdx)));
    end
    x = x(idx);  g = g(idx);  y = y(idx);

    %{
    % Estimate AR(4) for x_t to extract structural shock \varepsion_{1t} without constant
    X_lagged = lagmatrix(x, 1:4); % Create lag matrix
    X_lagged = X_lagged(5:end, :);
    x_regress = x(5:end);
    beta_x = (X_lagged' * X_lagged) \ (X_lagged' * x_regress); % OLS
    eps1 = x_regress - X_lagged * beta_x; % Extract ε_{1t}, 500
    eps1_full        = NaN(size(x));     % same length as W
    eps1_full(5:end) = eps1;             % first 4 rows are NaN, 504
    %}


    % Construct W_t (x_t, g_t, y_t)
    W = [x g y];
    
    % Correctly assign state indicator for each dgp_type
    switch dgp_type
        case 'TVAR' 
        S = [0; (y(1:end-1) > 1)]; % based on y_{t-1}

        case 'MSVAR'
        S = logistic_trans(W,0.5);
        % test 
        % fprintf('Recession frequency: %.1f%%\n', 100*mean(S==0));
        % fprintf('Expansion frequency: %.1f%%\n', 100*mean(S==1));

        case 'STVAR'
        % more extreme logistic functions to distict from TVAR
        kappa_val = 5; % or 2 or 10, experiment
        center_val = median(y,"omitmissing");  % or try mean(y), 1, 0.5, etc.
        S = smooth_trans(W, kappa_val, center_val);
    end
    
    % Find indices for each regime
    R_indices = find(S == 0); % Recession (S_{t-1} = 0)
    E_indices = find(S == 1); % Expansion (S_{t-1} = 1)

    % Separate regimes using find()
    W_R = W(R_indices, :); % Recession regime data
    W_E = W(E_indices, :); % Expansion regime data

    %{
    % If want to use VECM with nonstationary data, refer to below
    % But it violates state-dependent environment due to trends
    % Estimate reduced form VAR(4) parameters
    if data_trans == "level"
        % Estimate VECM 
        [k_R, A_R, C_R, Sigma_R] = estimateVECM(W_R, lag_type);
        [k_E, A_E, C_E, Sigma_E] = estimateVECM(W_E, lag_type);
  
    else
        [k_R, A_R, Sigma_R] = estimateVAR(W_R);
        [k_E, A_E, Sigma_E] = estimateVAR(W_E);
    % Estimate structural coefficient matrices C_R and C_E
    % lower–triangular version: C_lower = estimateC(W_state, eps1_state, true);
    C_R = estimateC(W_R, Sigma_R);
    C_E = estimateC(W_E, Sigma_E);
    %  Adjust sign in C
    end
    %}

    % Estimate reduced form VAR(4) parameters
    [k_R, A_R, Sigma_R] = estimateVAR(W_R);
    [k_E, A_E, Sigma_E] = estimateVAR(W_E);
    % Estimate structural coefficient matrices C_R and C_E
    % lower–triangular version: C_lower = estimateC(W_state, eps1_state, true);
    C_R = estimateC(W_R, Sigma_R);
    C_E = estimateC(W_E, Sigma_E);

    % To save gamma 
    gamma = calibrate_gamma(W, 4, 24, 0.0001, 0.45);

    % Store parameters
    params.k_RF_R = k_R; params.k_RF_E = k_E; 
    params.k_R = C_R*k_R; params.k_E = C_E*k_E;
    params.A_R = A_R; params.A_E = A_E;  % [] for VECM, coefficients for VAR
    params.C_R = C_R; params.C_E = C_E;
    %{
    % ----------------------------------------
    % test purpose: exaggerated difference
    params.C_R(2:end,:) = params.C_R(2:end,:) * 10;
    params.C_E(2:end,:) = params.C_E(2:end,:) * (-10); 
    %----------------------------------------
    %}
    params.Sigma_R = Sigma_R; params.Sigma_E = Sigma_E; 
    params.gamma = gamma;
end