function Y_b = bootstrap_VAR_sim(Bc,By,Y,Res_b,nlags)
% Generate bootstrap sample Y_b given residual draws
% Inputs:
%   Bc     = constant term (n×1)
%   By     = lag coefficient matrices (n×n×p)
%   Y      = original data (T×n)
%   Res_b  = bootstrap residuals (T-nlags × n)
%   nlags  = VAR order
%
% Output:
%   Y_b    = bootstrap sample (T×n)

    [T,n] = size(Y);
    Y_b = zeros(T,n);

    % initialize with actual history
    Y_b(1:nlags,:) = Y(1:nlags,:);

    % simulate forward
    for t = nlags+1:T
        yhat = Bc;
        for j = 1:nlags
            yhat = yhat + By(:,:,j) * Y_b(t-j,:)';
        end
        Y_b(t,:) = yhat' + Res_b(t-nlags,:);  % add bootstrap residual
    end
end
