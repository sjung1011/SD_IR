function Bh_init=setPriorMean_VAR(BVARcoeffs,h,modelSpec)
%initialize prior on powers of VAR posterior mean BVARcoeffs=[(n*nP)xn] (no
%constant)

%unpack basics
n  =modelSpec.modelSize; 
nL =modelSpec.nVARlags;
nP =modelSpec.nBLPlags;

%companion form
Ai=BVARcoeffs; %do not include constant

A=zeros(n*nL,n*nL); 
A(1:n,:)=Ai'; A(n+1:end,1:n*(nL-1))=eye(n*(nL-1));

%power of VAR coefficients
Bh =A^h; 
Bh =[zeros(1,n); Bh(1:n,:)']; %constant centered around zero


if nL > nP
    Bh_init=Bh(1:n*nP+1,:);
else
    Bh_init=[Bh; zeros(n*(nP-nL),n)];
end
