function [Bc,Br,Bx,By,Sigma,Sxx,w, y, X, x, Res,V] = LP_emp_contmp(Y,recurShock,respV,nlags,nhorizon)
% data for LP routine
nv = size(Y,2);
%nT = size(Y,1);
y  = Y(:,respV); % response variable
x  = lagmatrix(Y(:,recurShock), nhorizon); % impulse variable, contmp
% contemperaneous except for shock and lagged Y
if recurShock ==1
    w  = [ lagmatrix(Y(:,2:end), nhorizon) , lagmatrix( Y , (1:nlags) + nhorizon ) ];
else % recurShock == 8
    w  = [ lagmatrix(Y(:,1:(recurShock - 1)), nhorizon) , lagmatrix( Y , (1:nlags) + nhorizon ) ];
end
 
y = y((nlags + nhorizon + 1):end);
x = x((nlags + nhorizon + 1):end);
w = w((nlags + nhorizon + 1):end, :); 
X = [ones(size(x,1),1), x, w];

 % least-squares LP
[Beta,Sigma,Sxx,Res] = LS(y,X);

% store LP coefficients
Bc = Beta(1); % constant
Bx = Beta(2); % impulse variable
Br = Beta(3:nv+1); % contemperaneous control
By = Beta(nv+2:end); % lagged controls
By = reshape(By,[nv,nlags]);

% standard errors 
spec.se = 'white';

[T, N] = size(X);
XX_i = eye(N)/(X'*X); % inverse
if strcmp(spec.se, 'homoskedastic')
    V = sum(Res.^2)/(T-N)*XX_i;
elseif strcmp(spec.se, 'white')
    % White heteroskedasticity-consistent std errors
    Q = X'*diag(Res.^2)*X;
    V = XX_i*Q*XX_i * T/(T-N);
elseif strcmp(spec.se, 'nw')
    % Newey-West std errors
    M = floor(0.75*T^(1/3)); % Andrews (1991)
    % M = floor(1.3*sqrt(T)); % Lazarus, Lewis, Stock, Watson (2018)
    Q1 = zeros(N,N,T,M);
    for j = 1:M
        for i = j+1:T
            Q1(:,:,i,j) = (1-j/(M+1))*e(i)*e(i-j)...
                *(X(i,:)'*X(i-j,:) + X(i-j,:)'*X(i,:));
        end
    end
    Q = X'*diag(Res.^2)*X + sum(Q1,[3,4]);
    V = XX_i*Q*XX_i * T/(T-N);
end 
V=(V+V')/2; % make sure V is symmetric
V = V(recurShock,recurShock);
end