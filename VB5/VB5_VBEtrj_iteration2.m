function [mu,CovDiag0,CovDiag1,logDetLambda,Stt2,Sttau2,dmu2,dmx2]=...
    VB5_VBEtrj_iteration2(...
    datX,datT,datOne,datEnd,dim,EsPst,lambda_inv,alpha,EtrjOne,EtrjEnd,tau)
% [mu,CovDiag0,CovDiag1,logDetLambda,Stt2,Sttau2,dmu2,dmx2]=...
%    VB5_VBEtrj_iteration2(...
%    datX,datT,datOne,datEnd,dim,pst,lambda_inv,alpha,EtrjOne,EtrjEnd,tau)

% v 3: trying out Stfan Engblom's suggestions

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

%% assembly loop
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
  
%% inversion loop
if(0)
for nt=1:length(EtrjOne)
    Lind=EtrjOne(nt):EtrjEnd(nt);
    
    Ldiag0_nt=Ldiag0(Lind);
    Ldiag1_nt=Ldiag1(Lind);

    Lsp_nt=spdiags([ Ldiag1_nt Ldiag0_nt [0;Ldiag1_nt(1:end-1)]],-1:1,length(Lind),length(Lind));
    Sigma_nt=inv(Lsp_nt);
    CovDiag0_nt=(diag(Sigma_nt,0));
    CovDiag1_nt=[(diag(Sigma_nt,1)+diag(Sigma_nt,-1))/2;0];
    rMax=abs(diag(Lsp_nt));
    logDetLambda_nt=sum(log(rMax))+log(det(spdiags(1./rMax,0,length(Lind),length(Lind))*Lsp_nt));
    
    % explicit inversion and determinant computation     
    %[Cov0,Cov1,logDet]=triSym_triInv_rescale(Ldiag0_nt,Ldiag1_nt);          
    % compare: OK!
    %disp(num2str([max(abs(CovDiag0_nt-Cov0)./CovDiag0_nt) ...
    %              max(abs(CovDiag1_nt-Cov1)./CovDiag1_nt) ...
    %              abs(logDetLambda_nt-logDet) ]))%...
    %%%warning('double work in VB5_VBEtrj_iteration2, needs optimization')
    
    [CovDiag0_nt,CovDiag1_nt,logDetLambda_nt]=triSym_triInv_rescale(Ldiag0_nt,Ldiag1_nt);
    CovDiag0(Lind,1)=CovDiag0_nt;
    CovDiag1(Lind,1)=CovDiag1_nt;
    
    logDetLambda=logDetLambda+logDetLambda_nt;
end
end
% try sparse solver for the whole trajetory at once
LAM=spdiags([ Ldiag1 Ldiag0 [0;Ldiag1(1:end-1)]],-1:1,TtotHidden,TtotHidden);
%MU=LAM\nu;
%disp(num2str(max(max(abs(MU-mu))),3)) % OK!
mu=LAM\nu;

if(0) % stefan's suggestions, inversion did not work very well.
% seems to work, but not very fast
a=eig(LAM);
detLog1=sum(log(a));
% disp(num2str( (detLog1-logDetLambda)/logDetLambda))

% invert with trucated Newtos interations
% input A, assumed to be sparse (or else...)

% inversion by Newtons method with built-in diagonal truncation
nnewt=10;
n=TtotHidden;
ndiagseff=1;

if(~exist('iLAM','var'))
    iLAM=spdiags([ CovDiag1 CovDiag0 [0;CovDiag1(1:end-1)]],-1:1,TtotHidden,TtotHidden);
    %iLAM = sparse(1:n,1:n,1./diag(LAM)); % initial guess
end

for k = 1:nnewt
  iLAM = 2*iLAM-iLAM*LAM*iLAM;

  % truncate to ndiagseff diagonals
  [i,j,s] = find(iLAM);
  dkeep = find(abs(i-j) <= ndiagseff);

  % assemble next round
  iLAM = sparse(i(dkeep),j(dkeep),s(dkeep),n,n);
  % Detta fungerar då Newton är självkorrigerande. Pröva!
end

end


% inverting the whole matric in one go? nope, wrong answer
%[C0,C1,LDET]   =triSym_triInv_rescale_trjWise(Ldiag0,Ldiag1,EtrjOne,EtrjEnd,length(EtrjOne));
%[cC0,cC1,lLDET]=         triSym_d1Inv_trjWise(Ldiag0,Ldiag1,EtrjOne,EtrjEnd,length(EtrjOne));
%disp(num2str([max(abs(CovDiag0-C0)./CovDiag0) ...
%              max(abs(CovDiag1-C1)./CovDiag1) ...
%              abs((logDetLambda-LDET)/logDetLambda) ...
%              max(abs(CovDiag0-cC0)./CovDiag0) ...
%              max(abs(CovDiag1-cC1)./CovDiag1) ...
%              abs((logDetLambda-lLDET)/logDetLambda)]))
%%%% do it in one go with a custom algorithm:
[CovDiag0,CovDiag1,logDetLambda]=...
    triSym_d1Inv_trjWise(Ldiag0,Ldiag1,EtrjOne,EtrjEnd,length(EtrjOne));
%% precomputation loop
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

    







