function [irf,w, nT_GIC, nv_GIC, est_Var_HAC, est_Var_HAC_bc, ei_ls, ei_bc,SSR_ls,SSR_bc, TSS, SE, SE_bc, X_LPs] = IRF_LP(Y,recurShock,respV,nlags,nhorizons)
% Auxiliary function for estimating IRFs using least-squares LP

nT = size(Y,1);

% error checking
if (nhorizons + nlags) >= nT
    error("Number of horizons too large! No obs in sample!")
end

irf = zeros(1,nhorizons + 1);
nT_GIC = irf;
nv_GIC = irf;
bc =irf;
est_Var_HAC =irf;
est_Var_HAC_bc =irf;
SE=irf;
SE_bc=irf;







%----------------------------------------------------------------
% HAC estimator: LS LP & X struct for White SE (empI error bands)
%----------------------------------------------------------------
% settings 
ei_ls = struct();
ei_bc = struct();
SSR_ls = zeros(nhorizons + 1, 1);
SSR_bc = zeros(nhorizons + 1, 1);
TSS = zeros(nhorizons + 1, 1);
X_LPs = struct();

% go thru horizon 0 to horizon max
for h = 0:nhorizons     %nhorizons is 20, 
    field_name = sprintf('h%d', h);
    [~,~,Bx,~,~,~, w1, y, X, ~, Res] = LP(Y,recurShock,respV,nlags,h); % h-step ahead LP
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
    % for White SE (common for all LP-variants with the same lags)
    X_LPs.(field_name)= X ;
    % computing HAC for each horizon
    [est_Var_HAC(1,h+1),SE(1,h+1)] = hac_standard_error(Res,1,nT_GIC(1,h+1));
    if h==0
        w=w1; % Store control data vector for later bias correction if desired
    end
end
% Obtain 1 by h 'irf' and 'w' after LP. They are the inputs of LP_CorrectBias

%----------------------
% HAC estimator: BC LP
%----------------------
for h=0: nhorizons
    [IRF_corr, acf_corr] = LP_CorrectBias(irf, w); 
    T = size(w,1);

    if h>0
       bc(1,h+1) = (1/(T-h))*acf_corr(1:h)*IRF_corr(h:-1:1)';
    else 
       bc(1,h+1)=0; %h=0
    end
end

% Obtained 1 x 21 'bc'

for h=0: nhorizons
    field_name = sprintf('h%d', h); % Dynamically create the field name
    [~,~,~,~,~,~, ~, ~, ~, x, Res] = LP(Y,recurShock,respV,nlags,h);
    Res_bc = Res - x*bc(1,h+1); % only correct shock part
    ei_bc.(field_name) =Res_bc;
    Res_sq=Res_bc.^2;
    SSR_bc(h + 1, :) = sum(Res_sq, 1);
    [est_Var_HAC_bc(1,h+1),SE_bc(1,h+1)] = hac_standard_error(Res_bc,1,nT_GIC(1,h+1));
end

end