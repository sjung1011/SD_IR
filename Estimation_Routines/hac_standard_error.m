function [est_Var_HAC, SE_HAC] = hac_standard_error(sscore,KK,TT)

% sscore = X .* sscore if sandwitch estimator otherwise sscore is T by K residual vector.
% V_X_inv = (X'*X/T)^(-1);
% KK: dimention of dependent vars Y
% TT = size(Y,1) 
%% bandwidth: automatic data-dependent method

ss = sscore; % residual, T by 1 in LPs
nomin = 0; denomin = 0; % initial values
for kk = 1:KK 
    S_1 = ss(:,kk);
    rrho_hat1 = (S_1(1:end-1)' * S_1(1:end-1))^(-1) * S_1(1:end-1)' * S_1(2:end); 
    nomin = nomin + rrho_hat1^2/(1-rrho_hat1)^4; % scalar, updated from the initial value by adding to the accumulated previous.
    denomin = denomin + (1-rrho_hat1^2)^2/(1-rrho_hat1)^4; % scalar
end
ll = 1.8171*(nomin/denomin)^(1/3)*TT^(1/3); % scalar
ll = round(ll);
ll = min(round(TT),ll);
ll = max(2,ll); % ll is q in Kernel fcn

%% bandwidth:Newey-West
%{
ll = round(4*(TT/100)^(2/9));
%}

%% HAC variance
V_HAC = sscore' * sscore/TT; % Gamma_0, K by K
for hh = 1:ll-1 
    sscore_TT = sscore(hh+1:TT,:);
    sscore_TThh = sscore(1:TT-hh,:);
    GG_hat = sscore_TT' * sscore_TThh;
    V_HAC = V_HAC + 1/TT*kernel(hh,ll)*(GG_hat+GG_hat'); % hh=distance, ll=bandwidth
end

% est_Var_HAC = V_X_inv*V_HAC*V_X_inv/TT; if sandwitch
% SE_HAC = sqrt(diag(est_Var_HAC));

est_Var_HAC = V_HAC/TT; % As V_HAC added TT times in the loop, K by K 
SE_HAC = sqrt(est_Var_HAC); % Standard error
%% t-statistic
%t_stat = (bhat_ols-beta_null)./SE_HAC;