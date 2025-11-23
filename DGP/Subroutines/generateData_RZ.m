function [w, S] = generateData_RZ(T, params, p, T_burn, doBurn,dgp_type)
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
    % S      - T x 1 vector of regime indicators (0 or 1)
    

    T_sim = T + (doBurn * T_burn);      % total simulated length 1,000
    % Initialize variables
    w = zeros(T_sim, 3); % (x_t, g_t, y_t)
    S = zeros(T_sim, 1);

    % Generate initial values for the first p=4
    fixed_value =[0,0.19,1];
    w(1:p, :) = repmat(fixed_value, p, 1);
    %w(1:p, :) = randn(p, 3); % Initialize with random values

    % Simulate data
    for t = p+1:T_sim
        % Determine regime based on y_{t-1}
        switch dgp_type
            case 'TVAR'
                S(t) = (w(t-1, 3) > 1); % Expansion if y_t-1 > 1, else recession
            case 'MSVAR'
                gamma  = params.gamma;  
                zlag = w(t-1,3);                       % driver  y_{t-1}
                % row 1 (index 0) : recession (0)  → expansion (1)
                % row 2 (index 1) : expansion (1)  → recession (0)
                p01 = 1/(1 + exp(-(gamma(1,1) + gamma(1,2)*zlag)));  % 0→1
                p10 = 1/(1 + exp(-(gamma(2,1) + gamma(2,2)*zlag)));  % 1→0
                 if S(t-1)==0                           % currently in recession
                    S(t) = rand <  p01;                % flip to expansion?
                else                                   % currently in expansion
                    S(t) = rand >= p10;                % stay (1) or fall (0)
                end

            case 'STVAR'     
                kappa     = 5;                % slope of logistic weight
                c_thresh  = 0.9973; 
                g_t  = 1/(1 + exp(-kappa * (w(t-1,3) - c_thresh))); % logistic weight
                S(t) = (g_t > 0.5); 
        end

        % Select parameters for current regime
        if S(t) == 1
            k = params.k_RF_E; C = params.C_E; A = params.A_E; %Sigma = params.Sigma_E;
        else
            k = params.k_RF_R; C = params.C_R; A = params.A_R; %Sigma = params.Sigma_R;
        end

        % Construct lagged values for all 4 lags
        lagged_w = reshape(w(t-1:-1:t-p,:).', 3*p,1); % 3p*1 stacked as [w'{t-1}; … ; w'{t-4}]
        % Stack the A-matrices horizontally: 3x3x4 to 3 × 12 
        Ablock = reshape(A, 3, 3*p);          % same as cat(2,A(:,:,1),A(:,:,2),...,A(:,:,p))
        % Generate w_t  (row) 
        % C \ randn(3, 1) is equivalent to inv(C)*\ randn(3, 1) 
        % where randn(3,1)~ N(0,I) standard normal structural shock
        % generate w with reduced form equation
        w(t,:) = (k + Ablock * lagged_w + C \randn(3,1))';   %  1x3
        
        % U = chol(Sigma,'lower')*randn(3,1);
        % w(t,:) = (k + Ablock * lagged_w + U)';
    end
	
	
    % Drop burn-in if requested 
    if doBurn && T_burn>0 
        w = w(T_burn+1:end,:);
        S = S(T_burn+1:end);
    else % discard initial condition for CARs
        w = w(p+1:end, :);
        S = S(p+1:end);
    end
	
	% ----------- post-trim balance check ------------------------------------
	minObs  = 300;       % 100 observations per state
	
	counts = histcounts(S,[-0.5 0.5 1.5]);   % [ #recessions  #expansions ]
	if any(counts < minObs)
    error('generateData_RZ:TooFewObs', ...
          'After trimming only [%d %d] obs per regime; need ≥ %d.', ...
          counts(1), counts(2), minObs);
	end
	% --------------------------------------------------------------------
	
end
