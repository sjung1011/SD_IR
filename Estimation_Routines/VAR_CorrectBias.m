function [Bc, ByCorrect, Sigma, Sxx, Sigma_bc, Res_bc,SE_B] = VAR_CorrectBias(Y,nlags)
% Auxiliary function for estimating bias-corrected VAR coefficients

% first run least-squares VAR

nv = size(Y,2);
nT = size(Y,1);
[Bc,By,Sigma,Sxx] = VAR(Y,nlags);

% write VAR coefficients in companion form

VAR_companion_form = zeros(nv * nlags);
VAR_companion_form(1:nv,:) = reshape(By(:,:,:), [nv, nv*nlags]);
VAR_companion_form((nv+1):end, 1:(nv*(nlags-1))) = eye(nv*(nlags-1));

if max(abs(eig(VAR_companion_form))) > 1
    
    % if OLS is non-stationary then no correction (Kilian (1998), p. 220)
    
    ByCorrect = By;
    AyCorrect = reshape(ByCorrect, [nv, nv*nlags]);
else

    % if OLS is stationary, then standard bias correction in Pope (1989)
    
    Ay = reshape(By,[nv,nv * nlags]);
    CompAy = kron(diag(ones(1,nlags - 1),-1),eye(nv));
    CompAy(1:nv,:) = Ay; % VAR coefficients in companion form
    E = eig(CompAy);

    % compute several parts in the bias term
    
    sm = zeros(nv * nlags);
    for i = 1:(nv * nlags)
        sm = sm + E(i,1) * inv(eye(nv * nlags) - E(i,1) * CompAy');
    end

    Gamma0 = Sxx(2:end,2:end) - Sxx(2:end,1) * Sxx(2:end,1)';
    G = zeros(nv * nlags);
    G(1:nv,1:nv) = Sigma;

    b = G * (inv(eye(nv * nlags) - CompAy') + CompAy' / (eye(nv * nlags) - (CompAy') * (CompAy')) + sm) / Gamma0;
    Bias = - b / (nT - nlags);

    % correct the bias in VAR coefficients
    
    CompAyCorrect = CompAy - Bias;
    AyCorrect = CompAyCorrect(1:nv,:);
    ByCorrect = reshape(AyCorrect,[nv,nv,nlags]);

    % ensure stationarity of bias-corrected estimate

    indic_stat = 0;

    delta = 1;

    while indic_stat == 0

        % write bias-corrected VAR coefficients in companion form again
        
        VAR_companion_form = zeros(nv * nlags);
        VAR_companion_form(1:nv,:) = reshape(ByCorrect(:,:,:), [nv, nv*nlags]);
        VAR_companion_form((nv+1):end, 1:(nv*(nlags-1))) = eye(nv*(nlags-1));

        % check if stationary
        
        if max(abs(eig(VAR_companion_form))) < 1
            
            indic_stat = 1;
            
        else
            
            % ad-hoc adjustment as in Kilian (1998), p.220
            
            delta    = delta - 0.01;
            Bias_adj = delta * Bias;

            CompAyCorrect = CompAy - Bias_adj;
            AyCorrect = CompAyCorrect(1:nv,:);
            ByCorrect = reshape(AyCorrect,[nv,nv,nlags]);
            
        end
    end
end
    Beta = [Bc, AyCorrect]'; % [nv x 1, nv x (nv*p)]' = [nv x nk]', nk=(1+nv*p) 
    X = lagmatrix(Y,1:nlags); % regressors X
    Y = Y((nlags + 1):end,:); % nT x nv
    X = X((nlags + 1):end,:);
    X = [ones(size(X,1),1),X]; % nT x nk
    Res_bc = Y-X*Beta; % nT x nv
    
    % for BIC
    Sigma_bc = Res_bc' * Res_bc / (size(Y,1) - size(X, 2));
    % standard error for errorband in empirical
    Var_B =kron(inv(X'*X),Sigma_bc);
    SE_B =sqrt(diag(Var_B)); % (1+nv)*nv x 1
    SE_B =reshape(SE_B,size(X,2),nv); % (1+nv*p) x nv = 36 x 7
end