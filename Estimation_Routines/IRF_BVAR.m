function [irf_draws, BIC_bvar, Res, SE_B] = IRF_BVAR(VAR_coef_post_mean,VAR_coef_post_vce_inv,ShockVector,nhorizons,ndraw, Y)
% Auxiliary function for do posterior sampling of IRFs

nv = size(ShockVector,1);
nx = size(VAR_coef_post_mean,1) / nv;
nlags = (nx - 1) / nv;
irf_draws = NaN(nv,nhorizons + 1,ndraw);
VAR_coef_post_vce_sqrt = chol(inv(VAR_coef_post_vce_inv), 'lower');

% construct regressors X & placeholder for BIC
X = lagmatrix(Y,1:nlags);
Y = Y((nlags+1):end,:);
X = X((nlags+1):end,:);
X = [ones(size(X,1),1),X];
nT = size(Y,1);
Sigma_draws = NaN(nv,nv,ndraw);
Res_draws = NaN(nT,nv,ndraw);

% run multiple posterior draws
for idraw = 1:ndraw
    
    % make draws on VAR coef.
    this_VAR_coef = VAR_coef_post_mean + VAR_coef_post_vce_sqrt * randn(nx*nv, 1);
    
    % reshape to get By (reduced-form coef. on endogenous variables)
    this_Beta = reshape(this_VAR_coef, [nx, nv]);
    this_By = reshape(this_Beta(2:end,:),[nv,nlags,nv]); % lagged term
    this_By = permute(this_By,[3,1,2]);
    
    % compute IRF based on this draw
    irf_draws(:,:,idraw) = IRF_SVAR(this_By,ShockVector,nhorizons);
    
    % for BIC
    this_Res = Y-X*this_Beta;
    this_Sigma = this_Res'*this_Res/(nT-nx);
    Sigma_draws(:,:,idraw) = this_Sigma; 
    Res_draws(:,:,idraw) = this_Res; 
end
    
    % after draws, we use the information
    Sigma = mean(Sigma_draws,3);
    BIC =log(det(Sigma)) + (nx^2 * nlags + nx) * log(nT) / (nT);
    BIC_bvar = BIC*ones(1,nhorizons+1);
    Res =mean(Res_draws,3);
    
    % standard error for errorband in empirical
    Var_B =kron(inv(X'*X),Sigma);
    SE_B =sqrt(diag(Var_B)); % (1+nv)*nv x 1
    SE_B =reshape(SE_B,size(X,2),size(Y,2)); % (1+nv*p) x endoge vars = 36 x 7
end
