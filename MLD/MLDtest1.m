
clear

% generate some data
p0=1;
A=1;
pE=0;
D=1e6;
dt=10e-3;
tE=7e-3;

Ddt=D*dt;
lambda_tru=2*D*dt;
R  =tE/dt/6;
tau=tE/dt/2;
beta=tau*(1-tau)-R;

sigErr=25;

T=20;
numTrj=10;
dim=1;

[x,s,y]=VB5_diffusiveHMM_blur_detach(p0,A,pE,Ddt,R,tau,sigErr,dim,T,numTrj);
%x={x};y={y};

% add explicit error estimates
xv=cell(size(x));
for k=1:length(xv)
    xv{k}=[x{k} sigErr*ones(size(x{k}))];
end

dat=ML1_preprocess(xv,dim);
[logL,trj]=ML1_logLlambda(dat,lambda_tru,tau,R);

%% plot some stuff
figure(1)
clf
hold on
box on

errorbar(1:T,dat.x(dat.one(1):dat.end(1)),sigErr*ones(size(x{1})),'b')
plot(y{1},'r')
%errorbar(1:T+1,trj.mu(dat.Yone(1):dat.Yend(1)),sqrt(trj.CovDiag0(dat.Yone(1):dat.Yend(1))),'r')
plot(1:T+1,trj.mu(dat.Yone(1):dat.Yend(1)),'c','linew',1)
plot(1:T+1,trj.mu(dat.Yone(1):dat.Yend(1))+sqrt(trj.CovDiag0(dat.Yone(1):dat.Yend(1))),'c--')
plot(1:T+1,trj.mu(dat.Yone(1):dat.Yend(1))-sqrt(trj.CovDiag0(dat.Yone(1):dat.Yend(1))),'c--')

legend('data\pm err','true positions','est. positions\pm std')
%% estimate lambda

f=@(L)(-ML1_logLlambda(dat,L,tau,R));

[lambda_ML,fMin]=fminsearch(f,lambda_tru);

ll=lambda_tru*logspace(-0.3,0.3,30);
logL=ll;
for k=1:length(ll)
    logL(k)=-f(ll(k));
end

figure(2)
clf
hold on
plot(ll/lambda_tru,logL,'-k.')
plot(lambda_ML/lambda_tru,-fMin,'r*')
grid on

xlabel(' D / D_{true}')
ylabel('ln L')


