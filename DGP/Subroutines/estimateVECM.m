function [k, A, C, Sigma] = estimateVECM(W, numLags)
    % Estimate VECM parameters and structural C matrix for data W

    % Johansen cointegration rank test
    [h,~,~,~,~,~] = jcitest(W, 'lags', numLags-1, 'model','H1');
    r = find(h,1); 
    if isempty(r), r = 1; end

    % Estimate VECM
    vecmObj = vecm(W, 'NumLags', numLags-1, 'Rank', r);
    k = vecmObj.Constant;   % Intercept (if not estimated, zeros)
    Sigma = vecmObj.Covariance;
    C = chol(Sigma, 'lower');
    C = C ./ diag(C);       % Normalize diagonal to 1
    A = [];                 % Not defined in VECM

    % Optionally: return additional outputs as struct if needed
    % For simulation with VECM: you need k, C, Sigma, and the vecmObj
    % For simulation with VAR: use [k, A, C, Sigma] from estimateVAR
end
