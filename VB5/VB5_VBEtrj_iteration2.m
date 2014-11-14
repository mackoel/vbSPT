function [mu,CovDiag0,CovDiag1,logDetLambda,Stt2,Sttau2,dmu2,dmx2]=...
    VB5_VBEtrj_iteration2(...
    datX,datT,datOne,datEnd,dim,EsPst,lambda_inv,alpha,EtrjOne,EtrjEnd,tau)
% [mu,CovDiag0,CovDiag1,logDetLambda,Stt2,Sttau2,dmu2,dmx2]=...
%    VB5_VBEtrj_iteration2(...
%    datX,datT,datOne,datEnd,dim,pst,lambda_inv,alpha,EtrjOne,EtrjEnd,tau)

% internal variables
TtotHidden=sum(datT+1);   % total number of steps in hidden trajectory
TmaxTrjHidden=max(datT+1);   % length of the longest hidden trajectory
nu_nt=zeros(TmaxTrjHidden,dim);
%Ldiag0_nt=zeros(TmaxTrjHidden,1); % diagonals in the Lambda matrix
%Ldiag1_nt=zeros(TmaxTrjHidden,1); % off-diagonals in the Lambda matrix

% output variables
mu=zeros(size(datX,1),dim);
CovDiag0=zeros(TtotHidden,1);
CovDiag1=zeros(TtotHidden-1,1);
logDetLambda=0;

Stt2=zeros(size(datX,1),1);
Sttau2=zeros(size(datX,1),1);
dmu2=zeros(size(datX,1),dim);
dmx2=zeros(size(datX,1),dim);

% assembly loop
nu=zeros(TtotHidden,dim);
Ldiag0=zeros(TtotHidden,1);
Ldiag1=zeros(TtotHidden,1);

Mgamma=EsPst*lambda_inv';
Malpha=EsPst*alpha';
for nt=1:length(EtrjOne)
    X1=datOne(nt);
    XT=datEnd(nt);
    T=1+XT-X1; % length of current trajectory

    % parameter time averages in this trajectory
    %Mgamma_t=pst(X1:XT,:)*Epar.lambda_inv';
    %Malpha_t=pst(X1:XT,:)*Epar.alpha';
    
    nu_nt(1,:)=datX(X1,:)*Malpha(X1)*(1-tau);
    for k=1:dim
        nu_nt(2:T,k)=datX(X1+1:XT,k).*Malpha(X1+1:XT)*(1-tau)...
            +datX(X1:XT-1,k).*Malpha(X1:XT-1)*tau;
    end
    nu_nt(T+1,:)=datX(XT,:)*Malpha(XT)*tau;
    
    %Ldiag_nt=zeros(T+1,3);
    % diagonal elements
    
    Ldiag0_nt = [Mgamma(X1)+(1-tau)^2*Malpha(X1);
                 Mgamma(X1+1:XT)+Malpha(X1+1:XT)*(1-tau)^2+Mgamma(X1:XT-1)+Malpha(X1:XT-1)*tau^2;
                 Mgamma(XT)+Malpha(XT)*tau^2];

    % off-diagonal elements
    Ldiag1_nt    = [Malpha(X1:XT)*tau*(1-tau)-Mgamma(X1:XT);0];    
    
    Lind=EtrjOne(nt):EtrjEnd(nt);
    nu(Lind,:)=nu_nt(1:T+1,:);    
    Ldiag0(Lind,1)=Ldiag0_nt(1:T+1);
    Ldiag1(Lind,1)=Ldiag1_nt(1:T+1);
end
  
% inversion loop
for nt=1:length(EtrjOne)
    Lind=EtrjOne(nt):EtrjEnd(nt);
    
    Ldiag0_nt=Ldiag0(Lind);
    Ldiag1_nt=Ldiag1(Lind);
    nu_nt    =nu(Lind,:);
    
    Lambda_nt=zeros(T+1,T+1);
    for kk=1:T
        Lambda_nt(kk,kk+1)=Ldiag1_nt(kk);
    end
    Lambda_nt=Lambda_nt+Lambda_nt'+diag(Ldiag0_nt);
    
    % invert
    Sigma_nt=inv(Lambda_nt);    
    mu_nt=Lambda_nt\nu_nt(1:T+1,:);
    rMax=abs(diag(Lambda_nt));
    logDetLambda=logDetLambda+sum(log(rMax))+log(det(diag(1./rMax)*Lambda_nt));

    % put stuff back in place
    mu(Lind,:)=mu_nt;

    CovDiag0_nt=(diag(Sigma_nt,0));
    CovDiag1_nt=[(diag(Sigma_nt,1)+diag(Sigma_nt,-1))/2;0];    
    CovDiag0(Lind,1)=CovDiag0_nt;
    CovDiag1(Lind,1)=CovDiag1_nt;
    
    % explicit inversion and determinant computation     
    [Cov0,Cov1,logDet]=triSym_triInv_rescale(Ldiag0_nt,Ldiag1_nt);          
    warning('double work in VB5_VBEtrj_iteration2, needs optimization')
end

% precomputation loop
for nt=1:length(EtrjOne)
    Lind=EtrjOne(nt):EtrjEnd(nt);    
    Xind=datOne(nt):datEnd(nt);       % indices in measured positions

    mu_nt      =mu(Lind,:);
    CovDiag0_nt=CovDiag0(Lind);
    CovDiag1_nt=CovDiag1(Lind);
    
    T=length(Xind);
    dmu2(Xind,:)=diff(mu_nt,1,1).^2;
    dmx2(Xind,:)=(datX(Xind,:)-(1-tau)*mu_nt(1:T,:)-tau*mu_nt(2:T+1,:)).^2;
    Stt2(Xind,1)  =CovDiag0_nt(1:T,1)+CovDiag0_nt(2:T+1,1)-2*CovDiag1_nt(1:T,1);
    Sttau2(Xind,1)=(1-tau)^2*CovDiag0_nt(1:T,1)...
                      +tau^2*CovDiag0_nt(2:T+1,1)...
              +2*tau*(1-tau)*CovDiag1_nt(1:T,1);    
end

    







