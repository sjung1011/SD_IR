function [IR, theta, gamma, SSigma, nT_GIC, nv_GIC, ei_pen] = locproj(y, x, w, H_min, H_max, r, lambda)
% Function for estimating IRF and coefficients in penalized LP
    % Reference: Smooth/penalized local projection (Barnichon & Brownlees, 2019)
    
    %%% Input %%%
    % y:       response variable
    % x:       impulse variable
    % w:       controls (contemperaneous and lagged)
    % H_min:   minimum horizon
    % H_max:   maximum horizon
    % r:       order of finite difference operator
    % lambda:  penalty strength
    
    %%% Output %%%
    % IR:    impulse response
    % theta: penalized coef, for Xb. Warning: correspond to b_k in our paper
    % gamma: unpenalized coef, for W. Warning: correspond to \zeta_h and \phi_{h,l} in our paper
    
    % Design data matrix
    [B, Xb, W, Y_resw, Xb_resw] = locproj_design(y, x, w, H_min, H_max);
    
    % Compute impulse responses and coefficients using partitioned formula
    [IR, theta, gamma] = locproj_partitioned(y, B, Xb, W, Y_resw, Xb_resw, r, lambda);
    









    %-----------------------------
    % for GIC
    % for CVA: ei_pen, need to edit lenght of ei, which is residual.
    %-----------------------------
    % SSigma = zeros(1, H_max +1-H_min);
    nT_GIC = size(y,1);
    nv_GIC = size(theta,1) + size(gamma,1);
    ei_pen = struct();


    for h=1:H_max +1-H_min % H_max =20, H_min =0
        Res = y-Xb(:,:,h)*theta-W(:,:,h)*gamma(:,h);
        Res(isnan(Res)) = 0;
        SSigma(1,h) = Res'*Res/(nT_GIC-nv_GIC);
        
        if h==1 %h=0
            ei_pen.h0=Res; %T=196
        end

        if h==2
            ei_pen.h1=Res(h:end,:); %T=195
        end

        if h==3
            ei_pen.h2=Res(h:end,:);
        end

        if h==4
            ei_pen.h3=Res(h:end,:);
        end

        if h==5
            ei_pen.h4=Res(h:end,:);
        end

        if h==6
            ei_pen.h5=Res(h:end,:);
        end

        if h==7
            ei_pen.h6=Res(h:end,:);
        end

        if h==8
            ei_pen.h7=Res(h:end,:);
        end

        if h==9
            ei_pen.h8=Res(h:end,:);
        end

        if h==10
            ei_pen.h9=Res(h:end,:);
        end

        if h==11
            ei_pen.h10=Res(h:end,:);
        end

        if h==12
            ei_pen.h11=Res(h:end,:);
        end

        if h==13
            ei_pen.h12=Res(h:end,:);
        end

        if h==14
            ei_pen.h13=Res(h:end,:);
        end

        if h==15
            ei_pen.h14=Res(h:end,:);
        end

        if h==16
            ei_pen.h15=Res(h:end,:);
        end

        if h==17
            ei_pen.h16=Res(h:end,:);
        end

        if h==18
            ei_pen.h17=Res(h:end,:);
        end

        if h==19
            ei_pen.h18=Res(h:end,:);
        end

        if h==20
            ei_pen.h19=Res(h:end,:);
        end

        if h==21
            ei_pen.h20=Res(h:end,:);
        end

    end
end
