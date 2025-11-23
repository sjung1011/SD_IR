function irf = GIRF_from_VAR(By,Bc,Res,responseV,recurShock,H,delta)
% Compute deterministic GIRF given estimated VAR coefficients
%
% Inputs:
%   By, Bc   = VAR parameters
%   Res      = residuals (T-nlags × n), used for covariance
%   responseV, recurShock = indices
%   H        = horizon length - 1
%   delta    = shock size
%
% Output:
%   irf (H+1 × 1)

    [n,~,p] = size(By);

    % Companion form
    F_top = By(:,:,1);
    for j = 2:p
        F_top = [F_top By(:,:,j)];
    end
    F = [F_top; eye(n*(p-1)) zeros(n*(p-1),n)];
    K_big = [Bc ; zeros(n*(p-1),1)];
    sel   = eye(n, n*p);

    % Initial state = zeros
    st0 = zeros(n*p,1);

    % Structural impact: take reduced-form covariance
    Sigma = cov(Res);
    cholS = chol(Sigma,'lower');

    shockVec = cholS(:,recurShock) * delta;  % generalized impulse
    stS = st0; stB = st0;

    irf = zeros(H+1,1);
    for h = 0:H
        yB = sel*stB;
        yS = sel*stS;
        irf(h+1) = yS(responseV) - yB(responseV);

        % propagate
        stB = F*stB + K_big; % baseline
        stS = F*stS + K_big + [shockVec; zeros(n*(p-1),1)]*(h==0);
    end
end
