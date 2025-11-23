function [k, A, Sigma] = estimateVAR(W)
    % Estimate intercept (k) and VAR(4) coefficients (A)
    %
    % Inputs:
    % W - trimmed State-specific data (T_s x 3)
    %
    % Outputs:
    % k - Intercept vector (3 x 1): inv(C)*k in the equation
    % A - Lag coefficient matrix (3 x 3 x 4):  inv(C)*B(L) in the equation

    T_s = size(W, 1);
    X = lagmatrix(W, 1:4);
    X = [ones(size(X,1),1)  X]; % add constant
    X = X(5:end, :); % Remove NaN values
    Y = W(5:end, :);

    % Estimate coefficients
    [beta,Sigma,~,~] = LS(Y,X);
    %beta = (X' * X) \ (X' * Y); % =X\Y
    k = beta(1, :)'; % Intercept
    A = reshape(beta(2:end, :)', 3, 3, 4); % 3x3x4 lag matrix
    % Res = Y-X*beta; % residual
    % Sigma is cov of Res; % 3x3
end
