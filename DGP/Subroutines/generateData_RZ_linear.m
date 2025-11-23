function [w] = generateData_RZ_linear(T, params, p, T_burn, doBurn)
    % Generate synthetic data for state-dependent SVAR(4) using estimated parameters from data.
    %
    % Inputs:
    % T - Sample size, T_L for CARs; T_S to estimate them with LPs, VARs (after any burn-in)
    % params: state-dependent params
        % k = inv(C)*k, A=inv(C)*B(L) in the reduced-form equation   
    % p       : lag order
    % T_burn  : burn-in length (0 → none)
    % doBurn  : true---drop first T_burn rows
    %           false--- keep full path (set T_burn = 0 if unused)
    % Outputs:
    % w      - T x 3 matrix of simulated endogenous variables (x_t, g_t, y_t)
    
    

    T_sim = T + (doBurn * T_burn);      % total simulated length 1,000
    % Initialize variables
    w = zeros(T_sim, 3); % (x_t, g_t, y_t)
    
    % Generate initial values for the first p=4
    fixed_value =[0,0.19,1];
    w(1:p, :) = repmat(fixed_value, p, 1);
    %w(1:p, :) = randn(p, 3); % Initialize with random values

    % Simulate data
    for t = p+1:T_sim
        
        k = params.k_RF; C = params.C; A = params.A;
        % Construct lagged values for all 4 lags
        lagged_w = reshape(w(t-1:-1:t-p,:).', 3*p,1); % 3p*1 stacked as [w'_{t-1}; … ; w'_{t-4}]
        % Stack the A-matrices horizontally: 3x3x4 to 3 × 12 
        Ablock = reshape(A, 3, 3*p);          % same as cat(2,A(:,:,1),A(:,:,2),...,A(:,:,p))
        % Generate w with reduced form equation
        w(t,:) = (k + Ablock * lagged_w + C \randn(3,1))';   %  1x3
        
    end
    % Drop burn-in if requested 
    if doBurn && T_burn>0 
        w = w(T_burn+1:end,:);
        
    else % discard initial condition for CARs
        w = w(p+1:end, :);
        
    end
end