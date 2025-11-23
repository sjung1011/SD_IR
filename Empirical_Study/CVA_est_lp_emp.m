function [CVA_w] = CVA_est_lp_emp(Y, settings,ei_lp,ei_lp_contmp,recurShock)
% leave h out for stationary data.  

% preparations
IRF_hor = settings.est.IRF_hor;
nlags = settings.est.n_lag_short;


% all lps have the same lag orders 
% make struct array due to different dimensions in each horizon.
% xs, inv_w of each h are common for all modles as they have the same nlangs.
inv_w = struct(); 
xs = struct();

%------------------------
% LP  
%------------------------
% data setting
for h = 0:IRF_hor-1
x  = lagmatrix(Y(:,recurShock), h); % impulse variable, contmp
w  = lagmatrix( Y, (1:nlags) + h); 
% trimming 
x = x((nlags + h + 1):end);
w = w((nlags + h + 1):end, :); 
X = [ones(size(x,1),1), x, w];
field_name = sprintf('h%d', h); % Dynamically create the field name
xs.(field_name) = X; % T_h by nv_h
inv_w.(field_name) = inv(xs.(field_name)'*xs.(field_name)); % nv_h by nv_h
end
ee_ls =struct();

% recall residual structs and filling ee for each horizon
for h=0:IRF_hor-1
    field_name = sprintf('h%d', h);
    if h>0
        eye_mat =eye(2*h-1);
        T = size(ei_lp.(field_name),1); % 
        for t=h:1:T-h+1
            ehat_slec = ei_lp.(field_name)(t-h+1:t+h-1,:);
            xh = xs.(field_name)(t-h+1:t+h-1,:);
            Pth = xh*inv_w.(field_name)*xh';
            ehat_h_slec = inv(eye_mat-Pth)*ehat_slec;
            ee_ls.(field_name)(t-h+1,:) =ehat_h_slec(h,:) ;
        end
    else % h==0
        ee_ls.h0 = ei_lp.h0;
    end
end 

%------------------------
% LP_contmp  
%------------------------
% data setting
for h = 0:IRF_hor-1
    x  = lagmatrix(Y(:,recurShock), h); % impulse variable, contmp
    % contemperaneous except for shock and lagged Y
    if recurShock ==1
        w  = [ lagmatrix(Y(:,2:end), h) , lagmatrix( Y , (1:nlags) + h ) ];
    else % recurShock == 8
        w  = [ lagmatrix(Y(:,1:(recurShock - 1)), h) , lagmatrix( Y , (1:nlags) + h ) ];
    end
    x = x((nlags + h + 1):end);
    w = w((nlags + h + 1):end, :); 
    X = [ones(size(x,1),1), x, w];
    field_name = sprintf('h%d', h); % Dynamically create the field name
    xs.(field_name) = X; % T_h by nv_h
    inv_w.(field_name) = inv(xs.(field_name)'*xs.(field_name)); % nv_h by nv_h
end
ee_ls_contmp =struct();

% recall residual structs and filling ee for each horizon
for h=0:IRF_hor-1
    field_name = sprintf('h%d', h);
    if h>0
        eye_mat =eye(2*h-1);
        T = size(ei_lp_contmp.(field_name),1); % same as ei_bc
        for t=h:1:T-h+1
            ehat_slec = ei_lp_contmp.(field_name)(t-h+1:t+h-1,:);
            xh = xs.(field_name)(t-h+1:t+h-1,:);
            Pth = xh*inv_w.(field_name)*xh';
            ehat_h_slec = inv(eye_mat-Pth)*ehat_slec;
            ee_ls_contmp.(field_name)(t-h+1,:) =ehat_h_slec(h,:) ;
        end
    else % h == 0
        ee_ls_contmp.h0 = ei_lp_contmp.h0;
    end
end

%------------------------------------------------
% combine ee across models and obtain a1 (mxmxh)
%------------------------------------------------
ee_all = struct;
m = settings.est.n_weights_lp; % 2
a1 = zeros(m,m,IRF_hor); % quadratic term
a2 =zeros(m,1); % linear term (no linear penalty here)

for h = 0:IRF_hor-1
    field_name = sprintf('h%d', h);
    ee_all.(field_name) = [ee_ls.(field_name), ee_ls_contmp.(field_name)];
    a1(:,:,h+1) = ee_all.(field_name)'*ee_all.(field_name);
end

%------------------------------------
% Optimization, CVA_w (mxh)
%------------------------------------
% minimize 1/2*w'*Qw+f'*w (foc of Q^2*w)

CVA_w = zeros(m,IRF_hor);
for h=0:IRF_hor-1
      w0=ones(m,1)/m;
      options=optimset('LargeScale','off','Display','off');
      ww=quadprog(a1(:,:,h+1),a2,zeros(1,m),0,ones(1,m),1,zeros(m,1),ones(m,1),w0,options);
      ww=ww.*(ww>0); % Ensure non-negative weights
      CVA_w(:,h+1)=ww/sum(ww); % Normalize weights to sum to 1, (2 x H)
end


end