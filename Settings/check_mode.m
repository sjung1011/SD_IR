%% SPECIFIC SETTINGS FOR PRE-SPECIFIED ROBUSTNESS-CHECK MODES

% set up directory for robustness-check modes

mode_list   = {'nonlinear', 'linear','cumulative IRF','indicator','actual'};
save_mode_dir = mode_list{mode_type};

% rewrite some baseline settings in "shared.m" for different robustness check modes

switch mode_type

    case 1 % nonlinear
        
         % rewrite nothing and use all the settings in "shared.m"
   
    case 2 % linear

    case 3 % cumulative IRF
% -------------for empirical folder name--------------
    case 4 % indicator

    case 5 % actual
        
        % rewrite nothing and use all the settings in "shared.m"
        % cumulative IRF will be imputed in "run_combine.m"

end