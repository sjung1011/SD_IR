function [S, gg] = smooth_trans(W, kappa, c)

% ==============================================================
%  SMOOTH_TRANS  –  binary state S from a first-order logistic weight
% --------------------------------------------------------------
%  INPUTS
%    W       T×3 matrix  [x g y]
%    kappa:   positive scalar (slope of logistic; default 10)
%             Large kappa: steeper logistic → state series behaves almost like a threshold TVAR;
%             Small kappa: smoother weighting and more frequent flips around the boundary.
%    c:       scalar threshold  (default 1)
%
%  OUTPUT
%    S       T×1 binary vector  (0 = recession, 1 = expansion)
%    g       T×1 smooth weight g_t   (optional second output)
% --------------------------------------------------------------

    if nargin < 2 || isempty(kappa), kappa = 10; end   
    if nargin < 3 || isempty(c),     c     =  1; end   % matches TVAR

    y = W(:,3);                        % driver z_{t-1} = y_{t-1}
    T = size(y,1);

    % --- 1. compute logistic weight  g_t  -----------------------
    % Note: gg(1) is undefined because it needs y_{t-1}; pad with NaN.
    gg        = NaN(T,1);
    gg(2:T)   = 1 ./ ( 1 + exp( -kappa * ( y(1:end-1) - c ) ) ); % it should be S

    % --- 2. classify into binary state --------------------------
    S        = zeros(T,1);             % pre-allocate (recession by default)
    S(gg>0.5) = 1;                     % expansion when weight > 0.5
end
