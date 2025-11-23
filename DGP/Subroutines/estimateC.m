function C = estimateC(W, Sigma)
    % Estimate structural coefficient matrix C_t
    %
    % Inputs
    %     W              : T_s × 3 data matrix  [x_t , g_t , y_t](T_s x 3)
    %     epsilon1       : Structural shock ε_{1t} extracted from x_t's AR(4)
    %     useLowerTriang : (optional) logical flag
    %                    false  -> unrestricted C22  (default)
    %                    true   -> lower–triangular C22
    % Output:
    %     C : 3 × 3 structural matrix with ones on the diagonal
    
    % Take the Cholestky factor, then rescale rows so the diagonal =1
    L  = chol(Sigma , 'lower');   % L L' = Sigma ,   diag(L) > 0
    d  = diag(L);                 % positive scaling vector
    C  = L./d;                    % Row-rescaling puts ones on the diagonal, yielding C

    %{
    Columns 
    good = ~isnan(epsilon1);
    g = W(good, 2);
    y = W(good, 3);
    e1 = epsilon1(good);

    % Estimate C_{21,t} (2x1) using simultaneous regression
    C_21 = (e1 \ [g y])'; % returns 2×1 vector
    % C_21 = W(:,2:3) \ epsilon1; % g_t ~ ε_{1t}, and y_t ~ ε_{1t}
    
    % Estimate C_{22,t} (2x2) with diagonal restriction, but not lower triangle
    % Note that GHKP (2022) didn't give restrictions on the lower triangle matrix
    C_22 = eye(2); % Enforce diagonals as 1


    if useLowerTriang      % -------- lower–triangular restriction -------
        % estimate only the lower off-diagonal element (2,1)
        C_22(2,1) = (y' * g) / (y' * y);   %  g_t  ~  y_t
        % upper element (1,2) remains 0 by construction
    else                    % -------- unrestricted (symmetric example) --
        % estimate both off-diagonals freely
        C_22(1,2) = (g' * y) / (g' * g);   %  y_t  ~  g_t
        C_22(2,1) = (y' * g) / (y' * y);   %  g_t  ~  y_t
    end

    %% 3.  Assemble full 3 × 3 matrix
    % [ 1  0  0
    %  -C21(1)   C22(1,1)=1   C22(1,2)
    %  -C21(2)   C22(2,1)     C22(2,2)=1 ]
    C = [1,              0,              0;
        -C_21(1),        C_22(1,1),      C_22(1,2);
        -C_21(2),        C_22(2,1),      C_22(2,2)];
    %}
end
