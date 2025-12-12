function [CVA_w] = CVA_est_lp(Y, settings,ei_lp,ei_bclp,ei_blp,ei_pen,responseV,recurShock)
% leave h out for stationary data.  

% preparations
IRF_hor = settings.est.IRF_hor;
nlags = settings.est.n_lag;

% all lps have the same lag orders 
% make struct array due to different dimensions in each horizon.
% xs, inv_w of each h are common for all modles as they have the same nlangs.
inv_w = struct(); 
xs = struct();

for h=0:IRF_hor-1 % h=0:20 forecast horizon (nph in original coding)
[~,x,w] = LP_gen_data(Y,recurShock,responseV,nlags,h);
X = [ones(size(x,1),1), x, w];

field_name = sprintf('h%d', h); % Dynamically create the field name
xs.(field_name) = X; % T_h by nv_h
inv_w.(field_name) = inv(xs.(field_name)'*xs.(field_name)); % nv_h by nv_h
end


%------------------------
% LSLP and BCLP 
%------------------------
% ee_lp loop is checked. bc part is occuring some error. 
ee_ls =struct();
ee_bc =struct();

% recall residual structs and filling ee for each horizon
for h=0:IRF_hor-1
    field_name = sprintf('h%d', h);
    if h>0
        eye_mat =eye(2*h-1);
        T = size(ei_lp.(field_name),1); % same as ei_bc
        for t=h:1:T-h+1
            ehat_slec = ei_lp.(field_name)(t-h+1:t+h-1,:);
            xh = xs.(field_name)(t-h+1:t+h-1,:);
            Pth = xh*inv_w.(field_name)*xh';
            ehat_h_slec = inv(eye_mat-Pth)*ehat_slec;
            ee_ls.(field_name)(t-h+1,:) =ehat_h_slec(h,:) ;
        end
        for t=h:1:T-h+1
            ehat_slec = ei_bclp.(field_name)(t-h+1:t+h-1,:);
            xh = xs.(field_name)(t-h+1:t+h-1,:);
            Pth = xh*inv_w.(field_name)*xh';
            ehat_h_slec = inv(eye_mat-Pth)*ehat_slec;
            ee_bc.(field_name)(t-h+1,:) =ehat_h_slec(h,:) ;
        end
    elseif h == 0
        ee_ls.h0 = ei_lp.h0;
        ee_bc.h0 = ei_bclp.h0;
    end
end 



%------------------------
% BLP
%------------------------

% ee_BLP
ee_blp = struct(); 

for h=0:IRF_hor-1
    field_name = sprintf('h%d', h);
    if h>0
        eye_mat =eye(2*h-1);
        T = size(ei_blp.(field_name),1); % same as ei_bc
        for t=h:1:T-h+1
            ehat_slec = ei_blp.(field_name)(t-h+1:t+h-1,:);
            xh = xs.(field_name)(t-h+1:t+h-1,:);
            Pth = xh*inv_w.(field_name)*xh';
            ehat_h_slec = inv(eye_mat-Pth)*ehat_slec;
            ee_blp.(field_name)(t-h+1,:) =ehat_h_slec(h,:) ;
        end
    elseif h == 0
        ee_blp.h0 = ei_blp.h0;
    end
end 



%------------------------
% Pen-LP
%------------------------

ee_pen = struct(); 

for h=0:IRF_hor-1
    field_name = sprintf('h%d', h);
    if h>0
        eye_mat =eye(2*h-1);
        T = size(ei_pen.(field_name),1); % same as ei_bc
        for t=h:1:T-h+1
            ehat_slec = ei_pen.(field_name)(t-h+1:t+h-1,:);
            xh = xs.(field_name)(t-h+1:t+h-1,:);
            Pth = xh*inv_w.(field_name)*xh';
            ehat_h_slec = inv(eye_mat-Pth)*ehat_slec;
            ee_pen.(field_name)(t-h+1,:) =ehat_h_slec(h,:) ;
        end
    elseif h == 0
        ee_pen.h0 = ei_blp.h0;
    end
end 


%------------------------------------------------
% combine ee across models and obtain a1 (mxmxh)
%------------------------------------------------
ee_all = struct;
m = settings.est.n_weights_lp; %4
a1 = zeros(m,m,IRF_hor); % quadratic term
a2 =zeros(m,1); % linear term (no linear penalty here)

for h = 0:IRF_hor-1
    field_name = sprintf('h%d', h);
    ee_all.(field_name) = [ee_ls.(field_name), ee_bc.(field_name), ee_blp.(field_name), ee_pen.(field_name)];
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
     CVA_w(:,h+1)=ww/sum(ww); % Normalize weights to sum to 1, (4 x 21)
end

%else 
%    CVA_w =zeros(nmodels,IRF_hor);
%end

end