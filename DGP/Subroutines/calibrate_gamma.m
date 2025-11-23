function gamma = calibrate_gamma(W, Drec, Dexp, p01_hi, p10_lo)
%  W        :  T×3 matrix  [x g y]
%  Drec     :  average recession length   (e.g. 7  quarters)
%  Dexp     :  average expansion length   (e.g. 24 quarters)
%  p01_hi   :  Pr(expansion➜recession) when y is +1σ above its mean, e.g. 0.05
%  p10_lo   :  Pr(recession➜expansion) when y is −1σ below its mean, e.g. 0.03
%
%  RETURNS
%     gamma = [ α01  β01 ;     % expansion→recession row
%               α10  β10 ];    % recession→expansion row
%            (ready to plug into logistic_trans)

    y  = W(:,3);                      % driver  z_{t-1} = y_{t-1}
    mu = mean(y);
    sig= std(y);

    %--- target “baseline” switch probs (at z = mu) -----------------------
    p01_mu = 1/Dexp;                  % expansion to recession
    p10_mu = 1/Drec;                  % recession to expansion

    %--- helper ----------------------------------------------------------
    logit = @(p) log(p./(1-p));

    %--- row: expansion (1) to column: recession (0) ------------------------------
    gammaSlope_01 = (logit(p01_hi) - logit(p01_mu)) /  sig;
    gammaBar_01   =  logit(p01_mu) - gammaSlope_01*mu;

    %--- row: recession (0) to column: expansion (1) ------------------------------
    gammaSlope_10 = (logit(p10_lo) - logit(p10_mu)) / (-sig);
    gammaBar_10   =  logit(p10_mu) - gammaSlope_10*mu;

    gamma = [gammaBar_01  gammaSlope_01;
             gammaBar_10  gammaSlope_10];
end