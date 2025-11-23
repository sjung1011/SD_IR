function Yb = wild_bootstrap_VAR(Y, Res, p, T)
% ==================================================================
% Helper: wild (Rademacher) bootstrap using existing VAR residuals
% ==================================================================
% Y     : original sample (T × n)
% Res   : VAR residuals aligned with Y(p+1:T,:)
% p     : lag order
% T     : sample length (passed to avoid size(Y,1) each call)

    Yb = Y;                                  % keep first p obs unchanged

    % Rademacher multipliers ±1
    sgn = (randi([0 1], T-p, 1) * 2 - 1);

    % Bootstrap innovations
    Res_star = Res .* sgn;                   % (T-p) × n

    % Fitted values: Ydep - Res
    Yhat = Y(p+1:T,:) - Res;

    % Replace dependent block
    Yb(p+1:T,:) = Yhat + Res_star;
end