function [IRF, GIC, ei, SSR, SE_HAC] = BLP_est(Y, settings, responseV, normalizeV)
% State-dependent Bayesian LP for a single response variable, normalized by the shock
% Inputs:
%   Y:      T x n (n=3) data, already subsetted for regime
%   settings: structure
%   respV:  scalar, column of Y for response variable (e.g. 2 for g, 3 for y)
%   normalizeV:  scalar, column of Y for normalizing shock (e.g. 1 for x)
% Outputs:
%   IRF:    (IRF_hor x 1) normalized IRF for the response variable
%   extra_outputs: struct (can include diagnostics, GIC, SE, etc)
    % preparation
    IRF_hor = settings.est.IRF_hor;
    nlags = settings.est.n_lag;

    % --- Problem size
    [T, n] = size(Y);
    modelSpec.modelSize= n;
   
    % Setup model spec (as in your original BLP_est)
    modelSpec.nVARlags = nlags;
    modelSpec.nBLPlags = nlags;
    modelSpec.nLPlags  = nlags;
    modelSpec.nHorizons = IRF_hor - 1;
    modelSpec.bandsCoverage = 90;
    modelSpec.priorType = 'VAR';
    modelSpec.presample = true;
    modelSpec.identification = 'CHOL';
    modelSpec.shockSize = zeros(n,1);
    modelSpec.shockSize(normalizeV) = 1;
    modelSpec.shockVar = false(1,n);
    modelSpec.shockVar(normalizeV) = true;

    % Set up hyperpar options as in your code
    hyperPars.isrw      = true(1, n);
    hyperPars.lambda    = .4;
    hyperPars.lambdaC   = 1e5;
    hyperPars.lambdaP   = .4;
    hyperPars.miu       = 1;
    hyperPars.theta     = 2;
    hyperPars.alpha     = 2;
    hyperPriorsOptions.hyperpriors   = true;
    hyperPriorsOptions.Vc            = 1e5;
    hyperPriorsOptions.pos           = find(~hyperPars.isrw);
    hyperPriorsOptions.MNalpha       = [];
    hyperPriorsOptions.MNpsi         = false;
    hyperPriorsOptions.noc           = false;
    hyperPriorsOptions.sur           = false;
    hyperPriorsOptions.Fcast         = false;
    hyperPriorsOptions.hz            = modelSpec.nHorizons;
    hyperPriorsOptions.mcmc          = false;
    hyperPriorsOptions.Ndraws        = 1200;
    hyperPriorsOptions.Ndrawsdiscard = 200;
    hyperPriorsOptions.MCMCconst     = 1;
    hyperPriorsOptions.MCMCfcast     = false;
    hyperPriorsOptions.MCMCstorecoeff= false;
    hyperPriorsOptions.initialValues = hyperPars;
    GibbsOptions.iterations          = 1200;
    GibbsOptions.burnin              = 200;
    GibbsOptions.jump                = 1;
    hyperPriorsOptions.GibbsOptions  = GibbsOptions;

    % ------- Begin BLP estimation -------
    % (Copy your internal BLP_est code, replacing any usage of 'responseV' with respV, and 'normalizeV' with normalizeV)
    % ... [REMAINING CODE OMITTED HERE FOR BREVITY -- see your own code]
    % At the very end, instead of IRF = IRF_resp./IRF_normalize; IRF = IRF';
    
%unpack input structures
%model basics
nL = modelSpec.nVARlags; %lags in VAR prior
nP = modelSpec.nBLPlags; %lags in BLP
nH = modelSpec.nHorizons; 

shockS=modelSpec.shockSize;

%prior specification
priorType = modelSpec.priorType;

%identification scheme
iScheme=modelSpec.identification;

%variance of the VAR constant
lambdaC=hyperPriorsOptions.initialValues.lambdaC; %very large number

%find RW variables
isRandomWalk=hyperPriorsOptions.initialValues.isrw;

%Gibbs Sampler
GibbsOptions = hyperPriorsOptions.GibbsOptions;

nDraws =GibbsOptions.iterations; 
nBurn  =GibbsOptions.burnin; 
nJump  =GibbsOptions.jump;

%-VAR PRIORS---------------------------------------------------------------

%build matrix of relevant lagged Y
Ylag=NaN(T-nL,n*nL); %[Y_{t-1},...,Y_{t-p}]';
for j=1:nL
    
    Ylag(:,n*(j-1)+1:n*j)=Y(nL-j+1:end-j,:);
    
end

nT=size(Ylag,1);
W=Y(nL+1:end,:); %time adjusted dependend W_{t};
Y_init=mean(Y(1:nL,:)); %average of initial observations

%B are the projection coefficients @ different horizons (VAR @ h=0 & h=1)
%Sigma is the covariance of the projection residuals (VAR @ h=0 & h=1)
%Sigma~IW(S,a); vecB|Sigma~N(B,V); V=kron(Sigma,Omega); vecB=[n*(1+n*nLv)x1]
%set prior residual variance (Sigma) using univariate AR(1) residuals
sigmaj=NaN(n,1); 
for j=1:n
    
    sigmaj(j)=std( W(:,j)-[ones(nT,1) Ylag(:,j)]*...
        ([ones(nT,1) Ylag(:,j)]\W(:,j)) );
    
end

%IW prior for VAR residual variance
S_init=diag(sigmaj.^2);     %prior scale
a_init=n+2;                 %prior dof (E[Sigma_init]=S_init)

%Gaussian prior for VAR coefficients ~N(B_init,V_init) (equations in columns)
B_init=zeros(n*nL+1,n);
B_init(2:n+1,:)=diag(1.*isRandomWalk); %prior mean (coefficients)
nB=numel(B_init);

%projection set
YprojSet = [ones(nT,1) Ylag]; 

% * * * * * * * * * * * get optimal hyperparameters * * * * * * * * * * * %

parsAtMode=maxMLikelihoodVAR(W,YprojSet,B_init,sigmaj.^2,Y_init,hyperPriorsOptions);
lambda=parsAtMode.postmax.lambda; %overall tightness of NIW prior
parsAtMode.postmax.sigmahat  %%check size

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
%prior variance (coefficients)
Omega_init=inv(blkdiag(1/lambdaC,kron(diag(1:nL).^2,diag(sigmaj.^2))/lambda^2)); %think of it as inv(Xd'Xd); Xd dummy observations

%--VAR POSTERIOR-----------------------------------------------------------

%posterior variance (coefficients)
Omega_end = inv(inv(Omega_init)+YprojSet'*YprojSet); %Kadiyala&Karlsson(1997)

%posterior mean (coefficients)
B_end = parsAtMode.postmax.betahat;

%update parameters of the IW
%posterior scale
S_end=parsAtMode.postmax.sigmahat*(nT+a_init+n+1);

a_end=a_init+nT; %posterior degrees of freedom

%-------------------------------------
% settings for moving average, h=0,1
%-------------------------------------
nT_GIC = zeros(1,IRF_hor); % use Xh
nv_GIC = zeros(1,IRF_hor); % use XhprojSet
GIC = zeros(1,IRF_hor);


% when nlags =4, size(Uh_0)=196 * the number of w_t elements depending on experiments. 
Uh_0=W-YprojSet*B_end; % residual h=0, cannot get when h=1, thus I'll temporarily use it for both. 
Uh_0 = bsxfun(@minus, Uh_0,mean(Uh_0));


% for JMA case 2
% since residuals have different length due to h, 
% we should make a structure array to store variables. 
% Note that Operator '.*' is not supported for operands of type 'struct'
% Also note that ei.# is impossible variable name. 
ei =struct();

ei.h0 = Uh_0(:, responseV); %196x1
ei.h1 = Uh_0(2:end, responseV);

%contemporaneous transmission coefficients: 'CHOL'
Bzero=chol(cov(W-YprojSet*B_end)); 
Bzero=bsxfun(@rdivide,Bzero,diag(Bzero));
                
%IRF at posterior mean
irfs=nan(n,nH+1);
irfs(:,1)=shockS'*Bzero;                % h=0, on impact, 1 by n result row vector is used as the first column of irfs.
irfs(:,2)=shockS'*Bzero*B_end(2:n+1,:); %h=1

%if AR-based prior: save coefficients for initialization at future horizons
%use AR(1) coefficients to initialize priors on future horizons
BVARcoeffs=zeros(n*nL+1,n);
for j=1:n
    BVARcoeffs([1 j+1],j)=[ones(nT,1) Ylag(:,j)]\W(:,j);
end    

%----------------------------------------------------------%
%remove deterministic component
tempTrendVAR = handleVARtrend(W,Ylag,B_end(:),modelSpec);
x = tempTrendVAR.detrended;   %detrended data
%----------------------------------------------------------%

%store sampler's output
nRetainedDraws        =(nDraws-nBurn)/nJump; j=1;
StructuralIRFsCollect =NaN(n,nH+1,nRetainedDraws);
BzeroCollect          =nan(n,n,nRetainedDraws);

%shrinkage over horizons
optimalLambda=NaN(nH,1); optimalLambda(1)=lambda;

%sampling for VAR coefficients distribution
for i=1:nDraws
    
    Sigma_end=iwishrnd(S_end,a_end);    %draw from posterior IW
    
    %posterior for VAR coefficients ~N(B_end,V_end)    
    %variance
    V_end=kron(Sigma_end,Omega_end);
    
    %mean
    B_end=Omega_end*(Omega_init\B_init+YprojSet'*W);
    
    %check VAR stability
    isStableDraw=false;

    while ~isStableDraw
        
        vecB=B_end(:)+chol(V_end)'*randn(nB,1); %draw from posterior N
%         isStableDraw=checkVARstability(vecB,modelSpec);        
        isStableDraw=true;
        
    end
    
    %sample Bzero
            Bzero_=chol(cov(W-YprojSet*reshape(vecB,numel(vecB)/n,n) )); 
            Bzero_=bsxfun(@rdivide,Bzero_,diag(Bzero_));
   
    %compute and save IRFs
    if i>nBurn && mod(i,nJump)==0 %reduces dependence among draws
        
        %save Bzero for sampling @h>1
        BzeroCollect(:,:,j)=Bzero_;

        %on impact
        StructuralIRFsCollect(:,1,j)=shockS'*Bzero_; 
        
        %h=1
        vecB=reshape(vecB,numel(vecB)/n,n);
        StructuralIRFsCollect(:,2,j)=shockS'*Bzero_*vecB(2:n+1,:); j=j+1;
                                
    end
    
end

%-BLP AT HORIZON > 1-------------------------------------------------------

%loop over horizons
for h=2:nH   % nH=20
    
    
    if mod(h,10)==0 || h==nH

        system(['say horizon ' num2str(h)]);
    end

    %build relevant projection set (shift obs backward to match horizon)
    XhLag=NaN(T-(nP+h),n*nP);
    for j=h:nP+h-1
    
        XhLag(:,n*(j-h)+1:n*(j-h+1)) = x(nP+h-j+1:end-j,:);

    end  
        
    Xh=x(nP+1+h:end,:); nT=size(Xh,1);
    

    
    XhprojSet=[ones(nT,1) XhLag]; %XhprojSet=XhLag if no constant
    % for GIC
    nv_GIC(h+1) = size(XhprojSet,2); % h=2:20 not h=0,1
   
    Xh_init=mean(x(1:nP,:)); %mean of initial observations (should now be zero)
    
    
    %use univariate local projection to initialize scale (NW corrected)
    gammaU=NaN(n,1);
    
    for k=1:n
        
        projCoeffs=XhprojSet(:,[1 k+1:n:(n*nP+1)])\Xh(:,k);

        %univariate projection residuals
        u=Xh(:,k)-XhprojSet(:,[1 k+1:n:(n*nP+1)])*projCoeffs;

        u=bsxfun(@minus,u,mean(u)); GammaU=(u'*u)/nT; %HAC error (T*)covariance estimator
        
        nwLags=h-1; %u_{t+h|t} is MA(h-1)
        
        %HAC correction: prior scale
        nwWeights=(nwLags+1-(1:nwLags))./(nwLags+1);
        for j=1:nwLags

            gammaj=(u(j+1:nT,:)'*u(1:nT-j,:))/(nT-j);
            
            GammaU=GammaU+nwWeights(j)*(gammaj+gammaj');

        end

        gammaU(k)=sqrt(GammaU); %scalar
        
    end
    
    
    %prior on proj coeffs: mean
            %centered on relevant power of VAR coefficients
            Bh_init=setPriorMean_VAR(BVARcoeffs(2:end,:),h,modelSpec);
    nBh=numel(Bh_init);
    
    % * * * * * * * * * * get optimal hyperparameters * * * * * * * * * * %
    parsAtMode=maxMLikelihoodBLP(Xh,XhprojSet,Bh_init,gammaU.^2,Xh_init,hyperPriorsOptions,nP,h);
    lambdaP=parsAtMode.postmax.lambda;  %overall tightness of NIW prior
    optimalLambda(h)=lambdaP;

    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %
    
    %prior variance
    OmegaH_init=inv(blkdiag(1/lambdaC,kron(eye(nP),diag(gammaU.^2))/(lambdaP^2))); %think of it as inv(Xd'Xd); Xd dummy observations
        
    %posterior variance
    OmegaH_end=inv(inv(OmegaH_init)+XhprojSet'*XhprojSet); %Kadiyala&Karlsson(1997)

    %posterior mean
    B_end=parsAtMode.postmax.betahat;
    
    constant(h,:)=B_end(1,:);
    
    %IRF at posterior mean
    irfs(:,h+1)=shockS'*Bzero*B_end(2:n+1,:);

    %posterior scale
    Sh_end=parsAtMode.postmax.sigmahat*(nT+a_init+n+1);
    
    
    %correct for MA in proj residuals: sandwich covariance matrix for the
    %variance of the projection coefficients

    %-------------------------------------------------------
    %settings for moving average, projection residuals h>=2
    %-------------------------------------------------------
    Uh=Xh-XhprojSet*B_end;
    Uh=bsxfun(@minus,Uh,mean(Uh));  
    
    
    % GIC, CVA: residual set across horizonsm h>1
    if h==2
       ei.h2=Uh(:, responseV); % checked working, T=194
    end

    if h==3
       ei.h3=Uh(:, responseV); % checked working
    end

    if h==4
       ei.h4=Uh(:, responseV); % checked working 
    end
    
    if h==5
       ei.h5=Uh(:, responseV); % checked working 
    end

    if h==6
       ei.h6=Uh(:, responseV); % checked working
    end

    if h==7
       ei.h7=Uh(:, responseV); % checked working
    end

    if h==8
       ei.h8=Uh(:, responseV); % checked working
    end

    if h==9
       ei.h9=Uh(:, responseV); % checked working
    end

    if h==10
       ei.h10=Uh(:, responseV); % checked working
    end

    if h==11
       ei.h11=Uh(:, responseV); % checked working
    end

    if h==12
       ei.h12=Uh(:, responseV); % checked working
    end

    if h==13
       ei.h13=Uh(:, responseV); % checked working
    end

    if h==14
       ei.h14=Uh(:, responseV); % checked working
    end

    if h==15
       ei.h15=Uh(:, responseV); % checked working
    end

    if h==16
       ei.h16=Uh(:, responseV); % checked working
    end

    if h==17
       ei.h17=Uh(:, responseV); % checked working
    end

    if h==18
       ei.h18=Uh(:, responseV); % checked working
    end

    if h==19
       ei.h19=Uh(:, responseV); % checked working
    end

    if h==20
       ei.h20=Uh(:, responseV); % checked working
    end

    %
    nwWeights=(nwLags+1-(1:nwLags))./(nwLags+1);
    %HAC correction: posterior scale (mean of posterior IW distribution)
    for l=1:nwLags
        Gammal=(Uh(l+1:nT,:)'*Uh(1:nT-l,:))/(nT-l);
        Sh_end=Sh_end+nwWeights(l)*(Gammal+Gammal');
    end
    
    %coefficients' distribution
    j=1;
    for i=1:nRetainedDraws
        %draw for posterior variance (error)
        SigmaH_end=iwishrnd(Sh_end,a_end);
        %posterior variance (coefficients)
        GammaH_end=kron(SigmaH_end,OmegaH_end);
        %draw from posterior (coefficients)
        vecB=B_end(:)+chol(GammaH_end)'*randn(nBh,1); %draw from posterior N
        vecB=reshape(vecB,nBh/n,n);
        %sample Bzero
        Bzero_=BzeroCollect(:,:,j); 
        %irf
        StructuralIRFsCollect(:,h+1,j)=shockS'*Bzero_*vecB(2:n+1,:); 
        j=j+1;           
    end
end
           
StructuralIRFsCollect =sort(StructuralIRFsCollect,3);

%bands coverage
bandSize=modelSpec.bandsCoverage;
uBound=bandSize+(100-bandSize)/2; uBound=uBound/100;
lBound=(100-bandSize)/2;          lBound=lBound/100;

%--------------------------
% computing GIC
%--------------------------

nv_GIC(:,1:2) = nv_GIC(3); 
SSR = zeros(IRF_hor, 1); % 21 x 1
est_Var_HAC =zeros(1,IRF_hor);
SE_HAC = zeros(1,IRF_hor);

% go thru horizons
for h = 0:IRF_hor-1
    field_name = sprintf('h%d', h); % Dynamically create the field name
    Res_sq = ei.(field_name).^2;
    SSR(h + 1, :) = sum(Res_sq, 1);
    nT_GIC(h+1) = size(ei.(field_name), 1); % Access the field and get the size
    [est_Var_HAC(:,h+1), SE_HAC(:,h+1)] = hac_standard_error(ei.(field_name),1,nT_GIC(h+1));
    GIC(:, h+1) = log(det(est_Var_HAC(:,h+1))) + nv_GIC(h+1) * log(nT_GIC(h+1)) / nT_GIC(h+1); % det is a scalar 
end
GIC=GIC';

    %----------------------
    IRF_resp = irfs(responseV, :);      % IRFs for response variable
    IRF_normalize = irfs(normalizeV, 1); % IRF for normalization (impulse) at h=0
    IRF = IRF_resp ./ IRF_normalize; % Normalized IRF
    IRF = IRF';                     % Output as (IRF_hor x 1)
end
