function [irf_paths,irf, band_lo, band_hi, Res] = GIRF_VAR_est_emp(Y, settings,responseV,recurShock,delta_in)
% -------------------------------------------------------------
%  LS VAR(p) paths by drawing future innovations independently for baseline and shocked paths
%  which is not strictly Koop et al (1996)'s GIRF 
%  Y         (T × n)   sample generated in the current MC draw
%  nlags                VAR order  (p in your notes)
%  recurShock             column index of shocked variable
%  responseV              column index of response variable
%  H                    horizon of IRF  (integer, H ≥ 0)
%  B                    # future-path simulations (e.g. 1000)
%
%  irf      (H+1 × 1)   mean GIRF for responseV to a delta shock
%  Res      same as the LSVAR
% -------------------------------------------------------------
nlags = settings.est.n_lag;
H = settings.est.IRF_hor-1;
B = settings.est.n_girf_paths;  

% ---------- 1.  estimate the VAR(p), decide delta  ----------
[Bc,By,~,~,Res,~,~,~,~] = VAR(Y, nlags);   % By: (n × n × nlags)
n  = size(Y,2);                            % # variables
k  = Bc;                                   % constant vector, n x 1

if nargin < 5 || isempty(delta_in)
    switch lower(settings.est.delta_mode)         % choose in settings
        case 'unit'          % structural IRF, Δ = 1
            delta = 1;
        case 'onesd'         % 1 s.d. of reduced-form residual
            delta = std(Res(:,recurShock));
        case 'twosd'         % 2 s.d. (non-linearities check)
            delta = 2*std(Res(:,recurShock));
    end
else
    delta = delta_in;        % user supplied in the call
end

% ---------- 2.  build companion representation -------------------
F_top = By(:,:,1); % By is nxnxp
for j = 2:nlags
    F_top = [F_top  By(:,:,j)];   % n × (n*p)
end
F = [ F_top;           
      eye(n*(nlags-1))       zeros(n*(nlags-1), n) ];
K_big = [k ; zeros(n*(nlags-1),1)]; % constant in companion matrix
sel    = eye(n, n*nlags);          % selection matrix, pick current y_t

% ---------- 3.  pre-allocation -----------------------------------
diffPaths = zeros(H+1, B);                      % response variable only
T         = size(Y,1);
epsPool   = Res - mean(Res);                    % recentre residuals (zero mean condition)
% size(unique(epsPool,'rows'))
% ---------- 4.  simulation over B paths --------------------------
for b = 1:B
    % 4.1  draw random starting point (history of length nlags)
    t0   = randi([nlags+1 , T-H]);
    lagY = flipud(Y(t0-nlags:t0-1,:))';         % n × nlags
    stB  = lagY(:);                             % baseline state vector
    stS  = stB;                                 % shocked state vector
    % 4.2  draw future innovations
    T_res = size(Res,1);                 % T - nlags
    % draw future innovations independently for baseline and shocked paths
    idxB = randi(T_res, H+1, 1);
    idxS = randi(T_res, H+1, 1);
    epsB = epsPool(idxB, :);
    epsS = epsPool(idxS, :);
    % impose shock only on shocked path at h=0
    epsS(1, recurShock) = epsS(1, recurShock) + delta;
    %{
    idx   = randi(T_res, H+1, 1);        % draw (H+1) x 1 uniformly within residual matrix
    epsB     = epsPool(idx, :);                 % baseline ε_t
    epsS     = epsB;                            % set counterfactual eps
    epsS(1, recurShock) = epsS(1, recurShock) + delta;  %   delta = 1 shock
    %}

    % 4.3  propagate both paths with companion form 
    for h = 0:H
        yB = sel*stB; %+ epsB(h+1,:)';
        yS = sel*stS; %+ epsS(h+1,:)';
        diffPaths(h+1,b) = yS(responseV) - yB(responseV);

        stB = F*stB + K_big + [epsB(h+1,:)'; zeros(n*(nlags-1),1)];
        stS = F*stS + K_big + [epsS(h+1,:)'; zeros(n*(nlags-1),1)];
    end
end

% ---------- 5.  average across the B paths -----------------------
irf_paths = diffPaths;
irf = mean(diffPaths, 2);                        % (H+1)×1
band_lo = prctile(diffPaths,5,2);
band_hi = prctile(diffPaths,95,2);
end
