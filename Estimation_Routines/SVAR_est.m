function [IRF,BIC_var, Res, SSR, TSS, SE_B] ...
    = SVAR_est(Y_raw,settings,responseV,recurShock,normalizeV,bias_corrected)
% Function for estimating IRFs using least-squares VAR or bias-corrected VAR

% preparations
IRF_hor = settings.est.IRF_hor;
nlags = settings.est.n_lag;

% estimate VAR

if bias_corrected == 0
    [Bc,By,Sigma,Sxx,Res,Beta,~,~,SE_B] = VAR(Y_raw,nlags); % no bias correction
    Sigma_bc = Sigma;
    Res_sq = Res.^2;
    SSR = sum(Res_sq,1)'; % (1 x n)' sum residuals across t 
    SSR_mat = repmat(SSR,1,IRF_hor); % n x h, copy sumRes as VARs has horizon-specific Res
    SSR = SSR_mat(responseV,:); % 1 x h
    SSR = SSR'; % h x 1, stored by model
    % error band in empirical: asymptotic variance of VAR coeff
    SE_B = SE_B(1+recurShock,responseV); % 1 x 1
    SE_B = repmat(SE_B,1,IRF_hor); % 1 x h as horizon constant se

else
    [Bc,By,Sigma,~,Sigma_bc,Res_bc,SE_B] = VAR_CorrectBias(Y_raw,nlags); % with bias correction
    Res = Res_bc;
    Res_sq = Res.^2; % Square each element in Res
    SSR = sum(Res_sq,1)'; % (1 x n)' sum residuals across t 
    SSR_mat = repmat(SSR,1,IRF_hor); % n x h copy sumRes as VARs has horizon-specific Res
    SSR = SSR_mat(responseV,:); % 1 x h
    SSR = SSR'; % h x 1
    % error band in empirical: asymptotic variance of VAR coeff
    SE_B = SE_B(1+recurShock,responseV); % 1 x 1
    SE_B = repmat(SE_B,1,IRF_hor); % 1 x h
end
% horizon constant TSS for VARs
Y_TSS = Y_raw((nlags + 1):end,:);
Ybar = mean(Y_TSS,1); % 1 x n
DY = Y_TSS - repmat(Ybar, size(Y_TSS,1), 1); % T x n
DY_sq=DY.^2;
TSS = sum(DY_sq,1)'; % (1 x n)' sum DY across t 
TSS_mat = repmat(TSS,1,IRF_hor); % n x h, copy sumDY_t as VARs has horizon-constant sample.
TSS = TSS_mat(responseV,:); % 1 x h
TSS = TSS'; % h x 1, stored by model

G = chol(Sigma, 'lower'); % Warning: correspond to matrix C in our paper
ShockVector = G(:,recurShock);

% estimate IRF

IRF = IRF_SVAR(By,ShockVector,IRF_hor - 1); % IRF to one unit of shock
IRF = IRF(responseV,:) / IRF(normalizeV,1); % normalize by response of normalization variable
IRF = IRF';

%----------------------------------------
% BIC: Note that it should be outside 'if-end' to be output. 
%----------------------------------------
% BIC_var
BIC_X = lagmatrix(Y_raw,1:nlags);
BIC_X = BIC_X((nlags+1):end,:);
BIC_X = [ones(size(BIC_X,1),1),BIC_X];
BIC_nT =size(BIC_X,1);
BIC_nv = size(BIC_X,2);

if bias_corrected == 0
    BIC_var = log(det(Sigma)) + (BIC_nv^2 * nlags + BIC_nv) * log(BIC_nT) / (BIC_nT);
    BIC_var = BIC_var*ones(1,IRF_hor);
else
    BIC_var = log(det(Sigma_bc)) + (BIC_nv^2 * nlags + BIC_nv) * log(BIC_nT) / (BIC_nT);
    BIC_var = BIC_var*ones(1,IRF_hor);
end 

end