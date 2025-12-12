%% Common Setup for All Estimation Methods

% unpack settings
IRF_hor = settings.est.IRF_hor;
nlags = settings.est.n_lag;

% location of variables: data and linear
respV.g = settings.est.resp_vars(1); %2
respV.y = settings.est.resp_vars(2); %3
normV.g = settings.est.shock;
normV.y = settings.est.shock; % 1 in w ( shock has a unit impact on itself, and all other variables’ responses correspond to that same unit shock.)
recurShock = settings.est.shock; % location of impulse variable
normalizeV = recurShock;
% location of variables: state dependent LP equation