% PROXY SVAR with instrumental variable identification of IRFs
% built on Mertens and Ravn (2012) assumes only one shock of interest
%
% miranda 2014 smirandaagrippino@london.edu

function res=ProxySVARidentification(rfInnovations,policyVarPos,iVar)
% inputs
% rfInnovations =[NxT] matrix of reduced form innovations
% policyVarPos  =position of the policy indicator within system
% iVar          =[Tx1] vector of instrument for identification of policy shock
%
% output (structure)
% B    =[Nx1] vector of responses of VAR innovations to structural policy shock
% psi  =estimated correlation between instument and structural policy shock
% L    =reliability of instrument
% e    =[1xT] vector of realized policy shock
% stat =for the coefficients of the regression of the policy innovation on the instrument.

[N,T]=size(rfInnovations);

%place policy indicator first
iP =ismember(1:N,policyVarPos); 

rfInnov =[rfInnovations(iP,:);rfInnovations(~iP,:)]; 
rfInnov =rfInnov';


%coefficients of regression on instrument
betaIV =kron(eye(N),[ones(size(iVar)) iVar])\rfInnov(:); 
betaIV =reshape(betaIV,length(betaIV)/N,N)';

beta   =betaIV(:,2)./betaIV(1,2); %\beta_{21}\beta_{11}^{-1}


%t and F stats (regression on instrument of relevant innovation)
tempX  =[ones(size(iVar)) iVar]; 
tempU  =rfInnov(:,1)-tempX*betaIV(1,:)'; 

SigmaB =kron(diag(diag((tempU'*tempU)/T)),inv((tempX'*tempX)/T))/T;
t_Stat =betaIV(1,:)'./sqrt(diag(SigmaB));

tempY  =tempX*betaIV(1,:)'-mean(rfInnov(:,1)); k=length(betaIV(1,:)')-1;
F_Stat =((tempY'*tempY)/k)/((tempU'*tempU)/(T-k-1));


%identification
SigmaU =cov(rfInnov);
Gamma  =beta(2:end)*SigmaU(1,1)*beta(2:end)'-(SigmaU(2:end,1)*beta(2:end)'+...
    beta(2:end)*SigmaU(2:end,1)')+SigmaU(2:end,2:end);

B      =sqrt(SigmaU(1,1)-(SigmaU(2:end,1)-beta(2:end)*SigmaU(1,1))'/Gamma*...
    (SigmaU(2:end,1)-beta(2:end)*SigmaU(1,1))); %\beta_{11}
B      =B.*beta; %first column of B (u_{t}=B*e_{t})


%relevance
SigmaMU =cov([iVar,rfInnov]); 
phi     =SigmaMU(1,2:end)/B';


%realized shock sequence
tempBinv =(SigmaU(2:end,1)-beta(2:end)*B(1)^2)'/(SigmaU(2:end,2:end)-beta(2:end)*B(1)^2*beta(2:end)');
Binv     =[inv(B(1)-tempBinv*B(2:end)) -inv(B(1)-tempBinv*B(2:end))*tempBinv]; %first row of B^{-1}
e        =rfInnov*Binv';


%realized shock sequence (alternative)
tempX =[ones(size(iVar)) rfInnov]; 
e2    =tempX*(tempX\iVar);


%reliability
iM  =(iVar~=0); Gam =(sum(iM)/T)\phi;
eSq =e.^2; 
zSq =(iVar - Gam*e).^2;

L   =(Gam^2*sum( eSq(iM)) + sum(zSq(iM))) \ Gam^2*sum(eSq(iM));


%load output
res.B =NaN(N,1); res.B(iP,:)=B(1,1); res.B(~iP,1)=B(2:end,:);

res.phi   =phi; 
res.L     =L; 
res.e     =e; 
res.e2    =e2;
res.tstat =t_Stat(end); 
res.fstat =F_Stat;
res.Binv  =Binv';