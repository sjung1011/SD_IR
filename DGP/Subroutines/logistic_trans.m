function [S] = logistic_trans(W, targetShare)
% LOGISTIC_TRANS  Draw state series S_t from a TVTP-MSVAR while
%                 lowering p01_hi until the recession share ≤ targetShare.
%
% INPUT
%   W            T×3 matrix [x g y]
%   targetShare  desired upper bound on mean(S==0)  (default 0.50)
%
% OUTPUT
%   S            T×1   regime indicator  (0 = recession)
%   gamma        2×2   calibrated γ-matrix that generated S
%   share        realised recession frequency

    if nargin < 2,  targetShare = 0.50;  end

    % ----------------  fixed calibration targets ----------------
    Drec  = 4;      % average length of recession spells
    Dexp  = 24;     % average length of expansion spells
    p10lo = 0.45;   % recession→expansion when y is –1σ

    % ----------------  driver and storage ----------------------
    T  = size(W,1);
    y  = W(:,3);
    maxIter = 12;                 % safety stop

    % ------------  start with a modest p01_hi ------------------
    p01_hi = 0.0001;                % expansion→recession when y is +1σ

    for iter = 1:maxIter

        % 1) calibrate γ with the current p01_hi
        gamma = calibrate_gamma(W, Drec, Dexp, p01_hi, p10lo);

        % 2) simulate the path S_t
        S    = zeros(T,1);  S(1) = 1;        % start in expansion
        for t = 2:T
            zlag = y(t-1);
            p01 = 1/(1+exp(-(gamma(1,1)+gamma(1,2)*zlag))); % 0→1
            p10 = 1/(1+exp(-(gamma(2,1)+gamma(2,2)*zlag))); % 1→0
            if S(t-1)==0
                S(t) = rand < p01;
            else
                S(t) = rand >= p10;
            end
        end


        share = mean(S==0);      % recession frequency
      %{  
        fprintf('iter %2d : p01_hi = %.5f  →  share = %.3f\n',...
               iter, p01_hi, share);

        if share <= targetShare
            fprintf(' ✓  target reached (share ≤ %.2f)\n', targetShare);
            return
        end
        %}
        
        % 3) tighten only p01_hi (factor < 1 lowers it)
        p01_hi = max(1e-6, p01_hi * 0.6);
    end

    warning('Stopped at iter %d: recession share = %.3f still above target %.2f.',...
            maxIter, share, targetShare);
end

%{
function S = logistic_trans(W)
% Input 
%   W =[x g y]


    T  = size(W,1);      % total sample length  (incl. burn-in)
    S  = zeros(T,1);     % pre-allocate, start in state 0
    y = W(:,3);
    gamma = calibrate_gamma(W, 4, 24, 0.00001, 0.50); 

    for t = 2:T
        zlag = y(t-1);       % z_{t-1} = y_{t-1}
        
        % --- 1. transition probabilities at time t ---------------
        p01 = 1 / ( 1 + exp( -(gamma(1,1) + gamma(1,2)*zlag) ) );   % 0 → 1
        p10 = 1 / ( 1 + exp( -(gamma(2,1) + gamma(2,2)*zlag) ) );   % 1 → 0
        
        % --- 2. draw next state, rand is u[0,1] -----------------------------------
        if S(t-1)==0                       % currently in recession
            S(t) = rand < p01;             % switch to expansion 
        else                                % currently in expansion
            S(t) = rand >= p10;            % stay or switch back 
        end
    end
%}
