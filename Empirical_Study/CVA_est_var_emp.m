function [CVA_w] = CVA_est_var_emp(Y,settings,ei_var,ei_bvar,responseV)
% run in run_app_uncertainty.m
% leave-h-out cross validation averaging (CVA)
% all vars have the same lag orders: LSVAR, BVAR 
% residual vector: ei_m is (nT x nv)

% preparations
IRF_hor = settings.est.IRF_hor;
nlags = settings.est.n_lag;

% dimensions of xs, inv_w are the same over the horizon due to the same number of obs (nT) 
X = lagmatrix(Y,1:nlags); % regressors X
Y_adj = Y((nlags + 1):end,:); % nT x nv time adj
X = X((nlags + 1):end,:);
xs = [ones(size(X,1),1),X]; % nT x nk, same across horizons
inv_w = inv(xs'*xs);

%------------------------
% LSVAR, BVAR
%------------------------
% the standard VAR residuals are obtained for all time periods once the VAR coefficients are estimated.
% Yt = B0+B1Yt-1+...+BpYt-p+e
% as they extrapolate IRF, standard VAR residuals is proper for our perpose instead of horizon-specific residuals. 

% picking the residual corresponding to responseV to compare the LPs.
ei_var = ei_var(:, responseV);
ei_bvar = ei_bvar(:, responseV);

% ei_m, xs, inv_w do not vary across horizons, while ee_m does. 
ee_ls = struct();
ee_bvar = struct();

for h=0:IRF_hor-1
    T = size(ei_var,1); % same for ei_bc and ei_bvar
    field_name = sprintf('h%d', h);
    if h>0
        eye_mat =eye(2*h-1);
        for t=h:1:T-h+1
            % Adjust residuals for LS
            ehat_slec = ei_var(t-h+1:t+h-1,:);
            xh = xs(t-h+1:t+h-1,:);
            Pth = xh*inv_w*xh';
            ehat_h_slec = inv(eye_mat-Pth)*ehat_slec;
            ee_ls.(field_name)(t-h+1,:) = ehat_h_slec(h,:);
            % Adjust residuals for BVAR
            ehat_slec = ei_bvar(t-h+1:t+h-1,:);
            ehat_h_slec = inv(eye_mat-Pth)*ehat_slec;
            ee_bvar.(field_name)(t-h+1,:) = ehat_h_slec(h,:);
        end
    else
        % Directly assign residuals for h = 0
        ee_ls.h0 = ei_var;
        ee_bvar.h0 = ei_bvar;
    end
end  


%------------------------------------
% combine ee and obtain a1 (mxmxh)
%------------------------------------

m = settings.est.n_weights_var; % number of weights for VAR approaches
a1 = zeros(m,m,IRF_hor); % quadratic term for weights should be m by m
a2 = zeros(m,1); % linear term (no linear penalty here) []
ee_all = struct;

for h = 0:IRF_hor-1
    field_name = sprintf('h%d', h);
    ee_all.(field_name) = [ee_ls.(field_name), ee_bvar.(field_name)]; % T by m
    a1(:,:,h+1) = ee_all.(field_name)'*ee_all.(field_name);
end

%------------------------------------
% Optimization,CVA_w (mxh)
%------------------------------------

CVA_w = zeros(m,IRF_hor);
for h=0:IRF_hor-1
      w0=ones(m,1)/m;
      options=optimset('LargeScale','off','Display','off');
      ww=quadprog(a1(:,:,h+1),a2,zeros(1,m),0,ones(1,m),1,zeros(m,1),ones(m,1),w0,options);
      ww=ww.*(ww>0); % Ensure non-negative weights
     CVA_w(:,h+1)=ww/sum(ww); % Normalize weights to sum to 1, (m x H)
end

%else 
%    CVA_w =zeros(nmodels_var,IRF_hor);
%end 

end