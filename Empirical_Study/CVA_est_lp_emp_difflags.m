function [CVA_w] = CVA_est_lp_emp_difflags(Y, settings, ei_lp, ei_lp_contmp, recurShock)
%   Cross-Validated Averaging weights for two Local-Projection variants
%   (plain LP and LP with contemporaneous controls).
%
%   Y            :  T × n data matrix
%   settings     :  struct with fields
%                     .est.IRF_hor
%                     .est.n_lag_lp             (lag length for LP)
%                     .est.n_lag_lp_contmp      (lag length for LP_contmp)
%                     .est.n_weights_lp         (m = 2 here)
%   ei_lp        :  struct of residual matrices (horizons h0 … hH-1)
%   ei_lp_contmp :  same for LP_contmp
%   recurShock   :  column index of the shock variable
%
%   Returns
%       CVA_w     :  m × H matrix of weights   (m = 2 models)

% ---------------------------------------------------------------------
% 0.  Settings
% ---------------------------------------------------------------------
settings.est.n_lag_lp         = settings.est.n_lag;   % plain LP lags
settings.est.n_lag_lp_contmp  = settings.est.n_lag_short;    % LP-contmp lags

H  = settings.est.IRF_hor;

% allow backward compatibility if new fields are absent
if isfield(settings.est,'n_lag_lp')
    p_lp = settings.est.n_lag_lp;
else
    p_lp = settings.est.n_lag;           % old single lag value
end

if isfield(settings.est,'n_lag_lp_contmp')
    p_ct = settings.est.n_lag_lp_contmp;
else
    p_ct = settings.est.n_lag;           % fallback
end

m = settings.est.n_weights_lp;           % = 2

% containers for design matrices and their inverses
xs_lp     = struct();   invw_lp     = struct();
xs_ct     = struct();   invw_ct     = struct();

% ---------------------------------------------------------------------
% 1.  Build X and X'X⁻¹ for each horizon, for **each** model
% ---------------------------------------------------------------------
Tfull = size(Y,1);

for h = 0:H-1
    fname = sprintf('h%d',h);

    % ------------- plain LP (lags = p_lp) ----------------------------
    x  = lagmatrix(Y(:,recurShock),h);
    wL = lagmatrix(Y,(1:p_lp)+h);
    rows = (p_lp + h + 1) : Tfull;       % keep obs with full lags
    X_lp = [ones(numel(rows),1), x(rows), wL(rows,:)];
    xs_lp.(fname)   = X_lp;
    invw_lp.(fname) = pinv(X_lp'*X_lp);

    % ------------- LP-contmp  (lags = p_ct) --------------------------
    x  = lagmatrix(Y(:,recurShock),h);

    if recurShock == 1
        w0 = lagmatrix(Y(:,2:end),h);            % drop the shock col
    else
        w0 = lagmatrix(Y(:,1:(recurShock-1)),h); % vars before the shock
    end

    wL = lagmatrix(Y,(1:p_ct)+h);
    rows = (p_ct + h + 1) : Tfull;
    X_ct = [ones(numel(rows),1), x(rows), w0(rows,:), wL(rows,:)];
    xs_ct.(fname)   = X_ct;
    invw_ct.(fname) = pinv(X_ct'*X_ct);
end

% ---------------------------------------------------------------------
% 2.  Project residual stacks →  ễ( t , h )
% ---------------------------------------------------------------------
ee_lp = struct();        ee_ct = struct();

for h = 0:H-1
    fname = sprintf('h%d',h);

    if h == 0
        ee_lp.h0 = ei_lp.h0;
        ee_ct.h0 = ei_lp_contmp.h0;
        continue
    end

    eyeH = eye(2*h-1);

    % ---------- plain LP ----------
    T_e = size(ei_lp.(fname),1);
    k   = size(ei_lp.(fname),2);
    ee_lp.(fname) = NaN(T_e-2*h+1,k);

    for t = h : T_e-h
        e_seg = ei_lp.(fname)(t-h+1 : t+h-1,:);
        xh    = xs_lp.(fname)(t-h+1 : t+h-1,:);
        P     = xh * invw_lp.(fname) * xh';
        proj  = (eyeH - P) \ e_seg;      %  (2h-1) × k
        ee_lp.(fname)(t-h+1, :) = proj(h, :);   % middle row only
    end

    % ---------- LP-contmp ----------
    T_e = size(ei_lp_contmp.(fname),1);
    ee_ct.(fname) = NaN(T_e-2*h+1,k);

    for t = h : T_e-h
        e_seg = ei_lp_contmp.(fname)(t-h+1 : t+h-1,:);
        xh    = xs_ct.(fname)(t-h+1 : t+h-1,:);
        P     = xh * invw_ct.(fname) * xh';
        proj  = (eyeH - P) \ e_seg;
        ee_ct.(fname)(t-h+1, :) = proj(h, :);
    end
end

% ---------------------------------------------------------------------
% 3.  Assemble quadratic forms   a1(:,:,h) = ễ' ễ
% ---------------------------------------------------------------------
a1 = zeros(m,m,H);
a2 = zeros(m,1);    % no linear penalty

for h = 0:H-1
    fname = sprintf('h%d',h);

    % equalise sample sizes
    T1 = size(ee_lp.(fname),1);
    T2 = size(ee_ct.(fname),1);
    T  = min(T1,T2);

    E = [ee_lp.(fname)(1:T,:), ee_ct.(fname)(1:T,:)];   % T × 2

    a1(:,:,h+1) = E' * E;
end

% ---------------------------------------------------------------------
% 4.  Solve the quadratic program  (non-neg., sum=1)
% ---------------------------------------------------------------------
opts = optimset('LargeScale','off','Display','off');
CVA_w = zeros(m,H);

for h = 1:H
    w0 = ones(m,1)/m;
    w  = quadprog(a1(:,:,h), a2, zeros(1,m), 0, ones(1,m), 1, ...
                  zeros(m,1), ones(m,1), w0, opts);
    w  = max(w,0);           % numeric safeguard
    CVA_w(:,h) = w / sum(w); % normalise
end
end
