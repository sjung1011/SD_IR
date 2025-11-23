function params = estimateParams_RZ_linear(filename, sample_size, dgp_type)
    % Estimate parameters for state-dependent SVAR(4) model from real data
    %
    % Inputs:
    % filename     : path to RZDAT.xlsx
    % sample_size  : "long"  (1889q1–2015q4)  or "short" (1947q1–2015q4)
    % dgp_type     : 'LSVAR'
    % Outputs:
    % params - SVAR (C) and reduced-form VAR parameters (A,k)

    % ---------- Load data and construct series
    sheetname = 'rzdat';
    TBL = readtable(filename,'Sheet',sheetname,'PreserveVariableNames',true);
    % 
    pgdpLag  = lagmatrix(TBL.pgdp,1);
    ynormLag = lagmatrix(TBL.rgdp_pott6,1);
    x = TBL.news ./ (ynormLag .* pgdpLag);                     % military news
    g = TBL.ngov ./ (TBL.pgdp .* TBL.rgdp_pott6);              % government spending
    y = TBL.rgdp ./ TBL.rgdp_pott6;                       % real GDP
    
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
    
    % Not assign state indicator for each dgp_type
    
    % Estimate reduced form VAR(4) parameters
    [k, A, Sigma] = estimateVAR(W);
    
    % Estimate structural coefficient matrices C
    C = estimateC(W, Sigma);

    % Store parameters
    params.k_RF = k;  
    params.k = C*k; 
    params.A = A; 
    params.C = C; 
    params.Sigma = Sigma; 
end