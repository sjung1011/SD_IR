function [irf,irf_draws,nT_GIC, nv_GIC, est_Var_HAC, ...
    ei_ls, SSR_ls, TSS] ...
= IRF_LP_emp_contmp(Y,recurShock,respV,nlags,nhorizons)
% Auxiliary function for estimating IRFs using least-squares LP

nT = size(Y,1);

% error checking
if (nhorizons + nlags) >= nT
    error("Number of horizons too large! No obs in sample!")
end

irf = zeros(1,nhorizons + 1);
nT_GIC = irf;
nv_GIC = irf;
est_Var_HAC =irf;
% errorbands, direct normal draws
n_draw = 5e3;
irf_draws = zeros(nhorizons + 1, n_draw);
% settings 
ei_ls = struct();
SSR_ls = zeros(nhorizons + 1, 1);
TSS = zeros(nhorizons + 1, 1);


% go thru horizon 0 to horizon max
for h = 0:nhorizons     %nhorizons is 20, 
    field_name = sprintf('h%d', h);
    [~,~,Bx,~,~,~, ~, y, X, ~, Res,V] = LP_emp_contmp(Y,recurShock,respV,nlags,h); % h-step ahead LP
    irf(1,h + 1) = Bx;
    nT_GIC(1,h + 1) = size (y,1);
    nv_GIC(1,h + 1) = size (X,2);
    ei_ls.(field_name) =Res;
    Res_sq=Res.^2;
    SSR_ls(h + 1, :) = sum(Res_sq, 1); % SSR is a scalar at each h 
    ybar = mean(y,1); % scalar
    Dy = y - repmat(ybar, size(y,1), 1); % T x 1
    Dy_sq=Dy.^2;
    TSS(h + 1, :)=sum(Dy_sq,1); 
    irf_draws(h+1,:) = normrnd(Bx,sqrt(V),n_draw,1);
    % computing HAC for each horizon, duplication in fact
    [est_Var_HAC(1,h+1),~] = hac_standard_error(Res,1,nT_GIC(1,h+1));
   
end

end