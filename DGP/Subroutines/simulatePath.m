function Y_path = simulatePath(Y_hist4, params, delta, regime)
    % Simulate the path of Y_{t+h} for h = 0,...,H
    %
    % Inputs:
    % Y_t_initial : 4×3  [ y_{t-4}; y_{t-3}; y_{t-2}; y_{t-1}]
    % params     - Struct containing model parameters
    % delta      - Size of the structural shock u_t
    % regime     - Initial regime (0 or 1)
    %
    % Output:
    % Y_path     - (H+1) x 3 matrix of simulated impulse response paths

    % Extract parameters for the specified regime
    if regime == 0
        k = params.k_RF_R; C = params.C_R; A = params.A_R;
    else
        k = params.k_RF_E; C = params.C_E; A = params.A_E;
    end

    % Simulation settings
    H = 20; % Impulse response horizon

    % Initialize storage for simulated paths
    Y_path = zeros(H+1, 3);
    
    % -------- h = 0  : generate y_t  --------------
    lagMat = Y_hist4(4:-1:1 , :);      % now row1=y_{t-1}, row4=y_{t-4}
    y0 = k;
    
    for L = 1:4
        y0 = y0 + A(:,:,L) * lagMat(L,:)';
    end
    
    u0        = randn(3,1);   
    u0(1)     = u0(1) + delta;    % add structural shock at h=0
    y0        = y0 + (C \ u0);
    
    Y_path(1,:) = y0.';       % store as row, h==0
      

   % -------- forward recursion h = 1,…,H -----------------------
    for h = 1:H
        if h < 4
            recent = Y_path(h:-1:1 , :);
            need   = 4 - size(recent,1);
            past   = Y_hist4(4:-1:4-need+1 , :); % fill with history
            lagMat = [recent ; past];        % 4×3
        else
            lagMat = Y_path(h:-1:h-3 , :);   % h ≥ 4
        end

        y_next = k;
        
        for L = 1:4
            y_next = y_next + A(:,:,L) * lagMat(L,:)';
        end
        u_next  = randn(3,1);                 % no extra δ after h = 0
        y_next  = y_next + (C \ u_next);
        
        Y_path(h+1,:) = y_next.';   % save row h
    end
end
