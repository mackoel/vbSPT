function [W,C,F]=VB5_VBEMiterator(W,dat,varargin)
%% [W,C,F]=VB5_VBEMiterator(W,dat,varargin)
%
% Perform VBEM iterations, on the VB structure W, with data (structure)
% X, until convergence. This version accepts d-dimensional data (T by d
% matrices), and uses a 'short forward' model, i.e., hidden states are
% associated with steps, not positions, and the hidden state at time t is
% associated with the particle motion on the interval t -> t+dt.
%
% options:
% 'estimate'   : if given, add a field est2, which contains memory and
%                computer intensive estimates, such as the Viterbi path.
% 'slim'       : if given, remove some potentially bulky fields to decrease
%                storage footprint of the model.
%
% 'maxIter',n  : run at most n VBE iterations.
% 'relTolF',tf : convergence criterion for relative change in likelihood bound.
% 'tolPar' ,tp : convergence criterion for M-step parameters.
%                iterate until
%                |dF/F| <= tf and max|dPar(i)/Par(i)| <= tp,
%                where Par(i) are all parameters for the variational
%                distributions. Default values are
%                maxIter=1000, relTolF=1e-8, tolPar=1e-2.
% 'outputLevel', {0,1,2}
%                0: no output, 1: display convergence, not progress (default),
%                2: display convergence measures for each iteration
%
% This function uses the mex file HMM_multiForwardBackward.mexXXX from
% HMMcore for the computer intensive nner loops, where XXX is a platform
% dependent extension. Please refer to HMMcore/compile_code.m to
% (re)compile for your system. 

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB5_VBEMiterator, variational EM iterations in the vbSPT package
% =========================================================================
% 
% Copyright (C) 2014 Martin Lindén
% 
% E-mail: bmelinden@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by linking or combining it
%  with Matlab or any Matlab toolbox, the licensors of this Program grant you 
%  additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.

%% start of actual code
% process options
C.maxIter=1000;
C.relTolF=1e-8;
C.tolPar=1e-2;
do_estimates=false;
do_slim=false;
displayProgress=false;
displayExit=true;

if(nargin>2)        % then parse options
    k=1;            % argument counter
    kmax=nargin-2;  % stopping criterion
    while(k<=kmax)
        option=varargin{k};
        if(strcmpi(option,'estimate'))
            do_estimates=true;
            k=k+1;
        elseif(strcmpi(option,'slim'))
            do_slim=true;            
            k=k+1;
        elseif(strcmpi(option,'maxIter'))
            if(~isempty(varargin{k+1}))
                C.maxIter=varargin{k+1};
                if(~isnumeric(C.maxIter) || C.maxIter~=round(C.maxIter) || C.maxIter<=0)
                    error('VB1_VBEMiter: maxIter option must be followed by a positive integer.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'relTolF'))
            if(~isempty(varargin{k+1}))
                C.relTolF=varargin{k+1};
                if(~isnumeric(C.relTolF) || C.relTolF<=0)
                    error('VB1_VBEMiter: relTolF option must be followed by a positive number.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'tolPar'))
            if(~isempty(varargin{k+1}))
                C.tolPar=varargin{k+1};
                if(~isnumeric(C.tolPar) || C.tolPar<=0)
                    error('VB1_VBEMiter: tolPar option must be followed by a positive number.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'outputLevel'))
            if(~isempty(varargin{k+1}))
                outputLevel=varargin{k+1};
                if(~isnumeric(outputLevel) || ~ismember(outputLevel,[0 1 2]))
                    error('VB1_VBEMiter: outpuLevel option must be followed by 0, 1, or 2.')
                end
                if(outputLevel==0) % adjust output settings
                    displayExit=false;
                elseif(outputLevel==2)
                    displayProgress=true;
                end
            end
            k=k+2;
        else
            error(['VB1_VBEMiter: option ' option ' not recognized.'])
        end
    end
end
%% remove old estimates
if(isfield(W,'est'))
    W=rmfield(W,'est');
end
if(isfield(W,'est2'))
    W=rmfield(W,'est2');
end
if(isfield(W,'Fterms'))
    W=rmfield(W,'Fterms');
end
%% construct convergence structure
C.iter    =0;
C.converged=false;
C.W0=W;
%% preprocess the data and insert default aggregation
W.T=dat.T;
W.dim=dat.dim;
dim=W.dim;
tau=W.param.blur_tau;
R  =W.param.blur_R;
beta=tau*(1-tau)-R;
if(beta~=0)
    warning('Non-zero beta not formally applicable to this algorithm. Numbers take on different meaning!')
end
if(beta<0)
   error('VB5_VBEMiterator: inconsistent blur parameters ( tau*(1-tau)-R<0 =')
end
%trjStart=[1 dat.end(1:end-1)+1];    % indices to start of trajectories
%trjEnd  =dat.end-1;                  % indices to last position in each trajectory,
% which is also the last hidden state interval
N=size(W.PM.wB,1);
W.N=N;
if(~isfield(W.M,'SA')) % add default state aggregation (no aggregation)
    W.M.SA=1:W.N;
end
%% initialize VBEM iterations
runMore=true;
C.exitStatus='';
Wm2=struct;Wm1=struct; %#ok
E0par=struct;
%% iterate
while(runMore)
    % keep a short history in case something goes wrong...
    Wm2=Wm1;Wm1=W; %#ok
    % try one round of VBEM iterations
    %% parameter M-step 1: W.E  + W.PM -> W.M, W.Epar
    if(isfield(W,'E')) % then actually do the update step        
        % transitions and initial conditions
        W.M.wPi = W.PM.wPi + W.E.wPi;
        wB=W.E.wA.*(1-eye(W.N)).*(W.PM.wB>0); % only allowed transitions included
        W.M.wa =  W.PM.wa  + [sum(wB,2) diag(W.E.wA)];
        W.M.wB =  W.PM.wB  + wB;        
        % emission model part, with aggregated states
        for a=1:max(W.M.SA)
            % all states in aggregate a gets emission statistics from
            % all states in the same aggregate
            ind=find(a==W.M.SA);
            W.M.nl(ind)  = W.PM.nl(ind)  + sum(W.E.nl(ind));
            W.M.cl(ind)  = W.PM.cl(ind)  + sum(W.E.cl(ind));            
            W.M.na(ind)  =                 sum(W.E.na(ind));
            W.M.ca(ind)  =                 sum(W.E.ca(ind));
        end
        % numerical evaluation of parameter expectation values, show
        % warning for potentiall badly conditioned parameters
        %%% really, only W.Epar.lambda_inv and W.Epar.alpha are used before
        %%% next M-step. Some time could be save omitting other integrals.
        [W.Epar.alpha,~,W.Epar.lambda_inv]=VB5_ag_vl_expectations_beta0(W);
        
        % check for problems
        isNanInf=(sum(~isfinite([W.M.wPi W.M.wa(1:end) W.M.wB(1:end) ...
                                 W.M.nl  W.M.cl W.M.na  W.M.ca ...
        W.Epar.alpha W.Epar.ln_alpha W.Epar.lambda_inv ...
        W.Epar.ln_lambda]))>1);
        if(isNanInf)
            error('VB5_VBEMiter:Mfield_not_finite','Nan/Inf generated in VBM step')
        end        
    end
    %% trajectory E-step : W.Es + W.M  -> W.Etrj
    if( isfield(W,'Es') && isfield(W,'M'))% && false)        
        % construct \nu
        Ttot=sum(dat.T+1);   % total number of steps in hidden trajectory
        nu=zeros(Ttot,dim);
        Ldiag=zeros(Ttot,3); % diagonals in the Lambda matrix
        for nt=1:length(W.Etrj.one)
            MU1=W.Etrj.one(nt);
            MUT=W.Etrj.end(nt)-1; % 1:T for hidden trajectory, which has T+1 points
            X1=dat.one(nt);
            XT=dat.end(nt);
            T=1+XT-X1; % length of current trajectory
            
            % parameter time averages in this trajectory
            Mgamma_t=W.Es.pst(X1:XT,:)*W.Epar.lambda_inv';
            Malpha_t=W.Es.pst(X1:XT,:)*W.Epar.alpha';
            
            for k=1:dim
                nu(MU1+1:MUT,k)=dat.x(X1+1:XT,k).*Malpha_t(2:T)*(1-tau)...
                    +dat.x(X1:XT-1,k).*Malpha_t(1:T-1)*tau;
            end
            nu(MU1,:)=dat.x(X1,:)*Malpha_t(1)*(1-tau);
            nu(MUT+1,:)=dat.x(XT,:)*Malpha_t(T)*tau;
            
            % diagonal elements
            Ldiag(MU1,2)      =Mgamma_t(1)+(1-tau)^2*Malpha_t(1);
            Ldiag(MU1+1:MUT,2)=Mgamma_t(2:T)+Malpha_t(2:T)*(1-tau)^2 ...
                +Mgamma_t(1:T-1)+Malpha_t(1:T-1)*tau^2;
            Ldiag(MUT+1,2)    = Mgamma_t(T)+Malpha_t(T)*tau^2;
            % off-diagonal elements
            Ldiag(MU1:MUT,1)    =Malpha_t*tau*(1-tau)-Mgamma_t;
            Ldiag(MU1+1:MUT+1,3)=Malpha_t*tau*(1-tau)-Mgamma_t;
        end
        % compute \Sigma(t,t), \Sigma(t,t+1)
        %Lambda=spdiags(Ldiag,-1:1,Ttot,Ttot);
        % invert Lambda block by block
        % Sigma1=inv(Lambda);
        % Sigma=spalloc(Ttot,Ttot,sum((W.Etrj.end-W.Etrj.one+1).^2));
        logDetLambda=0;
        
        mu=zeros(size(nu));
        CovDiag0=zeros(size(nu,1),1);
        CovDiag1=zeros(size(nu,1)-1,1);
        
        for nt=1:length(W.Etrj.one)
            tic
            Lind=W.Etrj.one(nt):W.Etrj.end(nt);
            
            Ldiag_nt=Ldiag(Lind,1:2);
            Lambda_nt=zeros(size(Ldiag_nt,1),size(Ldiag_nt,1));
            for kk=1:size(Ldiag_nt,1)-1
                Lambda_nt(kk,kk+1)=Ldiag_nt(kk,1);
            end
            Lambda_nt=Lambda_nt+Lambda_nt'+diag(Ldiag_nt(:,2));
            
            Sigma_nt=inv(Lambda_nt);
            nu_nt=nu(Lind,:);
            
            mu_nt=Lambda_nt\nu_nt;
            CovDiag0_nt=(diag(Sigma_nt,0));
            CovDiag1_nt=(diag(Sigma_nt,1)+diag(Sigma_nt,-1))/2;
            
            %Sigma(Lind,Lind)=Sigma_nt;
            
            mu(Lind,:)=mu_nt;
            CovDiag0(Lind,1)=CovDiag0_nt;
            CovDiag1(Lind(1:end-1),1)=CovDiag1_nt;
            
            
            rMax=abs(diag(Lambda_nt));
            %T_nt=length(Lind);
            logDetLambda=logDetLambda+sum(log(rMax))+log(det(diag(1./rMax)*Lambda_nt));
        end
        tic
        % determinant ln|Lambda| with scaling to
        % rMax=abs(diag(Lambda));
        % W.Etrj.logDetLambda=sum(log(rMax))+log(det(spdiags(1./rMax,0,Ttot,Ttot)*Lambda));
        W.Etrj.logDetLambda=logDetLambda;
        clear Lind rMax Lambda_nt logDetLambda %T_nt
        
        % collect averages for hidden trajectory
        %W.Etrj.mu=Lambda\nu;
        %W.Etrj.CovDiag0=full(diag(Sigma,0));
        %W.Etrj.CovDiag1=full(diag(Sigma,1)+diag(Sigma,-1))/2; % possible better to average of rounding errors?
        W.Etrj.mu=mu;
        W.Etrj.CovDiag0=CovDiag0;
        W.Etrj.CovDiag1=CovDiag1; % possible better to average of rounding errors?
        
        clear mu Covdiag1 CovDiag0 mu_nt Covdiag1_nt CovDiag0_nt Sigma_nt nu_nt Ldiag_nt
        clear nu Ldiag Lambda_nt MU1 MUT X1 XT T Ttot Mgamma_t Malpha_t
        clear Ldiag Sigma logDetLambda
    end    
    %% partial parameter M-step 2
    if(isfield(W,'Etrj') && isfield(W,'Es'))
        W.E.cl=zeros(1,N);
        W.E.ca=zeros(1,N);
        for nt=1:length(W.Etrj.one)
            MU1=W.Etrj.one(nt);
            MUT=W.Etrj.end(nt)-1; % 1:T for hidden trajectory, which has T+1 points
            X1=dat.one(nt);
            XT=dat.end(nt);
            T=1+XT-X1; % length of current trajectory
            
            X     = dat.x(X1:XT,:);
            MU    = W.Etrj.mu(MU1:MUT+1,:);
            Stt   = W.Etrj.CovDiag0(MU1:MUT+1);
            Sttp1 = W.Etrj.CovDiag1(MU1:MUT);
            W.E.cl= W.E.cl+(...
                dim/2*( Stt(1:T)+Stt(2:T+1)-2*Sttp1(1:T))...
                +1/2*sum(diff(MU,[],1).^2,2) ...
                )'*W.Es.pst(X1:XT,:);
            W.E.ca= W.E.ca+(...
                dim/2*((1-tau)^2*Stt(1:T)+tau^2*Stt(2:T+1)+2*tau*(1-tau)*Sttp1(1:T))...
                +1/2*sum((X-(1-tau)*MU(1:T,:)-tau*MU(2:T+1,:)).^2,2)...
                )'*W.Es.pst(X1:XT,:);
        end
        
        for a=1:max(W.M.SA)
            % all states in aggregate a gets emission statistics from
            % all states in the same aggregate. Only some expecation values
            % changed since last M-step
            ind=find(a==W.M.SA);
            W.M.cl(ind)  = W.PM.cl(ind)  + sum(W.E.cl(ind));            
            W.M.ca(ind)  =                 sum(W.E.ca(ind));
        end
        % check for problems
        isNanInf=(sum(~isfinite([W.M.cl W.M.ca]))>1);
        if(isNanInf)
            error('VB5_VBEMiter:Mfield_not_finite','Nan/Inf generated in VBM step 2')
        end
    end
    %% hidden state E-step, and complain/crash if it cannot be done
    % numerical evaluation of parameter expectation values, show
    % warning for potentiall badly conditioned parameters
    [W.Epar.alpha,W.Epar.ln_alpha,W.Epar.lambda_inv,W.Epar.ln_lambda,...
        W.Epar.KL_lv]=VB5_ag_vl_expectations_beta0(W);
    % check for problems
    isNanInf=(sum(~isfinite([W.Epar.alpha W.Epar.ln_alpha ...
        W.Epar.lambda_inv W.Epar.ln_lambda W.Epar.KL_lv]))>1);
        if(isNanInf)
            error('VB5_VBEMiter:Epar_not_finite','Nan/Inf generated before E step')
        end
    
    if( isfield(W,'M') && isfield(W,'Etrj'))
        % coupling matrix is the same for all trajectories
        % lnQ=psi(W.M.wA)-psi(sum(W.M.wA,2))*ones(1,W.N); % old version
        lnQ=zeros(W.N,W.N);
        wB0=sum(W.M.wB,2);
        wa0=sum(W.M.wa,2);
        for i=1:W.N
            lnQ(i,i)=psi(W.M.wa(i,2))-psi(wa0(i));
            for j=[1:(i-1) (i+1):W.N]
                lnQ(i,j)=psi(W.M.wa(i,1))-psi(wa0(i))...
                    +psi(W.M.wB(i,j))-psi(wB0(i));
            end
        end
        
        lnQmax=max(max(lnQ));
        Q=exp(lnQ-lnQmax);
        % for estimates of transition rates and dwell times: s(t)
        W.est.Q=Q;
        W.est.lnQ=lnQ;
        
        % assemble time-dependent weights
        lnH1=psi(W.M.wPi)-psi(sum(W.M.wPi)); % weights from initial state probability
        lnH0=dim/2*(-W.Epar.ln_lambda+W.Epar.ln_alpha);   % data-independent term, same for all t
        lnH =ones(size(dat.x,1),1)*lnH0;
        
        for nt=1:length(W.Etrj.one)
            MU1=W.Etrj.one(nt);
            MUT=W.Etrj.end(nt)-1; % 1:T for hidden trajectory, which has T+1 points
            X1=dat.one(nt);
            XT=dat.end(nt);
            T=1+XT-X1; % length of current trajectory
            MU=W.Etrj.mu(MU1:MUT+1,:);
            X=dat.x(X1:XT,:);
            Stt  =W.Etrj.CovDiag0(MU1:MUT+1);
            Sttp1=W.Etrj.CovDiag1(MU1:MUT);
            for j=1:N
                lnH(dat.one(nt),j)=lnH(X1,j)+lnH1(j); % start of each trajectory
                lnH(X1:XT,j)=lnH(X1:XT,j)... % dimension-independent data terms
                    -dim/2*W.Epar.lambda_inv(j)*(Stt(1:T)+Stt(2:T+1)-2*Sttp1)...
                    -dim/2*W.Epar.alpha(j)*((1-tau)^2*Stt(1:T)+tau^2*Stt(2:T+1)+2*tau*(1-tau)*Sttp1);
                for k=1:dim
                    lnH(X1:XT,j)=lnH(X1:XT,j)...
                        -1/2*W.Epar.lambda_inv(j)*(diff(MU(:,k))).^2 ...
                        -1/2*W.Epar.alpha(j)*(X(:,k)-(1-tau)*MU(1:T,k)-tau*MU(2:T+1,k)).^2;
                end
            end
        end
        lnHMax=max(lnH,[],2);
        H=zeros(size(lnH));
        for j=1:N % now we compute the actual pointwise emission contribution
            H(:,j)=exp(lnH(:,j)-lnHMax);
        end

        [lnZz,wA,W.Es.pst]=HMM_multiForwardBackward(Q,H,dat.end);
        % forward sweep normalization constant (same as VB3)
        lnZQ=(sum(dat.T-2))*lnQmax;
        lnZq=sum(lnHMax);
        
        % transition counts
        W.E.wA=wA;
        W.E.wPi=sum(W.Es.pst(dat.one,:),1);          % <s_1>=<\delta_{j,s_1}>_{q(s)}            
        
        % statistics for next parameter update
        W.E.nl= dim/2*sum(W.Es.pst,1);
        W.E.na=W.E.nl;

        W.E.cl=zeros(1,N);
        W.E.ca=zeros(1,N);
        for nt=1:length(W.Etrj.one)
            MU1=W.Etrj.one(nt);
            MUT=W.Etrj.end(nt)-1; % 1:T for hidden trajectory, which has T+1 points
            X1=dat.one(nt);
            XT=dat.end(nt);
            T=1+XT-X1; % length of current trajectory
           
            X     = dat.x(X1:XT,:);
            MU    = W.Etrj.mu(MU1:MUT+1,:);
            Stt   = W.Etrj.CovDiag0(MU1:MUT+1);
            Sttp1 = W.Etrj.CovDiag1(MU1:MUT);
            W.E.cl= W.E.cl+(...
                        dim/2*( Stt(1:T)+Stt(2:T+1)-2*Sttp1(1:T))...
                         +1/2*sum(diff(MU,[],1).^2,2) ...
                    )'*W.Es.pst(X1:XT,:);
            W.E.ca= W.E.ca+(...
                dim/2*((1-tau)^2*Stt(1:T)+tau^2*Stt(2:T+1)+2*tau*(1-tau)*Sttp1(1:T))...
                 +1/2*sum((X-(1-tau)*MU(1:T,:)-tau*MU(2:T+1,:)).^2,2)...
                )'*W.Es.pst(X1:XT,:);
        end
    else
        error('VB5_VBEMiterator: not enough model fields to perform the E-step.')
    end
        % check for problems
    isNanInf=(sum(~isfinite([W.E.nl W.E.cl W.E.na W.E.ca]))>1);
    if(isNanInf)
        error('VB5_VBEMiter:Efield_not_finite','Nan/Inf generated in VBEs step')
    end
    %% lower bound
    % hidden trajectory
    if(isfield(W.Etrj,'logDetLambda'))
        W.Fterms.trj=dim/2*sum(dat.T+1)*(1+log(2*pi))-dim/2*W.Etrj.logDetLambda;
        F=W.Fterms.trj;
    else
        F=0;
        W.Fterms.trj=0;
    end
    
    % hidden states
    F=F+lnZQ+lnZq+lnZz;
    W.Fterms.lnZQ=lnZQ;
    W.Fterms.lnZq=lnZq;
    W.Fterms.lnZz=lnZz;
    if(~isfinite(F))
        error('VB5_VBEM: F not finite (lnZ)')
    end    
    
    % KL divergence of transition probabilities of s(t).
    % Not changed from VB3 -> VB5
    KL_a=zeros(W.N,1);    
    if(W.N>1) % a is only defined if N>1
        wa0=sum(W.M.wa,2);
        ua0=sum(W.PM.wa,2);
        KL_a=gammaln(wa0)-gammaln(ua0)...
            -(wa0-ua0).*psi(wa0)-(...
            gammaln(W.M.wa(:,1))-gammaln(W.PM.wa(:,1))...
            -(W.M.wa(:,1)-W.PM.wa(:,1)).*psi(W.M.wa(:,1))...
            +gammaln(W.M.wa(:,2))-gammaln(W.PM.wa(:,2))...
            -(W.M.wa(:,2)-W.PM.wa(:,2)).*psi(W.M.wa(:,2)));
    end
    W.Fterms.aTerms=-KL_a;
    F=F-sum(KL_a);
    if(~isfinite(F))
        error('VB5_VBEM: F not finite (KL_a)')
    end
    clear wa0 ua0;    
    % jump probabilities
    KL_B=zeros(1,W.N);
    if(W.N>1) % B is only defined for N>1        
        for k=1:W.N
            %ind=setdiff(1:N,k); % only include non-diagonal elements
            ind=find(W.PM.wB(k,:)>0); % only include non-zero elements
            wB0=sum(W.M.wB(k,ind));
            uB0=sum(W.PM.wB(k,ind));
            KL_B(k)=gammaln(wB0)-gammaln(uB0)-(wB0-uB0)*psi(wB0)...
                +sum((W.M.wB(k,ind)-W.PM.wB(k,ind)).*psi(W.M.wB(k,ind))...
                -gammaln(W.M.wB(k,ind))+gammaln(W.PM.wB(k,ind)));
        end
    end
    W.Fterms.Bterms=-KL_B;
    F=F-sum(KL_B);
    if(~isfinite(F))
        error('VB5_VBEM: F not finite (KL_B)')
    end    
    clear wA0 uA0 ind;
    % KL divergence of initial state probability 
    u0Pi=sum(W.PM.wPi);
    w0Pi=sum(W.M.wPi);
    KL_pi=gammaln(w0Pi)-gammaln(u0Pi)...
          +sum((gammaln(W.PM.wPi)-gammaln(W.M.wPi))...
               +(W.M.wPi-W.PM.wPi).*(psi(W.M.wPi)-psi(w0Pi)));
    W.Fterms.piTerms=-KL_pi;
    F=F-KL_pi;
    if(~isfinite(F))
        error('VB3_VBEM: F not finite (KL_pi)')
    end
    
    % KL divergence of emission parameters
    KL_lv=W.Epar.KL_lv;
    %KL_gj= W.PM.ng.*log(W.M.cg./W.PM.cg)...
    %    -W.M.ng.*(1-W.PM.cg./W.M.cg)...
    %    -gammaln(W.M.ng)+gammaln(W.PM.ng)...
    %    +(W.M.ng-W.PM.ng).*psi(W.M.ng);
    %KL_aj= W.PM.na.*log(W.M.ca./W.PM.ca)...
    %    -W.M.na.*(1-W.PM.ca./W.M.ca)...
    %    -gammaln(W.M.na)+gammaln(W.PM.na)...
    %    +(W.M.na-W.PM.na).*psi(W.M.na);
    % remove duplicate terms in each aggregate
    for a=1:max(W.M.SA)
       ind=find(a==W.M.SA);
       KL_lv(ind(2:end))=0;
       %KL_gj(ind(2:end))=0;
       %KL_aj(ind(2:end))=0;
    end    
    W.Fterms.lvTerms=-KL_lv;
    %W.Fterms.alphaTerms=-KL_aj;
    F=F-sum(KL_lv);
    if(~isfinite(F))
        error('VB5_VBEM: F not finite (KL_lv)')
    end
    %% assembly of the free energy
    W.F=F;                
    if(~isfinite(W.F))
        error('VB5_VBEMiter:F_not_finite','Nan/Inf generated in lower bound')
    end
    %catch me
    %% catch potential errors
    %me.getReport
    %crashfile=sprintf('VB1_VBEMiterator_error_%f.mat',now);
    %save(crashfile)
    %runMore=false;
    %C.exitStatus=[C.exitStatus 'VB1_VBEMiterator generated error'];
    %error(['VB1_VBEMiterator: VB iterations returned error. Saved state to ' crashfile])
    %end
    C.iter=C.iter+1;
    %% do convergence check!
    if(isfield(Wm1,'F')) % then we can check convergence
        %% converge criterion in terms of relative changes in F and parameter values
        if(~isfinite(W.F)) % check for problem
            crashfile=sprintf('VB5_VBEMiterator_NaNInf_%f.mat',now);
            save(crashfile);
            error(['VB5_VBEMiterator found W.F=NaN or Inf. Saving state to ' crashfile])
        end
        C.dFrel=(W.F-Wm1.F)/abs(W.F); % convergence statistic for F
        % convergence statistics for parameters
        fM=fields(W.M);
        C.dPrel=-Inf;
        C.limitPar='';
        for k=1:length(fM)
            Pnew=W.M.(fM{k})(1:end);
            Pold=Wm1.M.(fM{k})(1:end);
            ind=find(Pnew~=0);
            dPrel=max(abs((Pnew(ind)-Pold(ind))./Pnew(ind)));
            if(dPrel>C.dPrel)
                C.limitPar=['M.' fM{k}];
                C.dPrel=dPrel;
            end
        end
        
        Fconverged=(abs(C.dFrel)<C.relTolF);
        Pconverged=C.dPrel<C.tolPar;
        if(Fconverged && Pconverged)
            C.exitStatus=['Converged normally after ' int2str(C.iter) ' iterations, relTolF = ' num2str(C.relTolF) ', tolPar = ' num2str(C.tolPar) ', (' C.limitPar ' limiting).'];
            C.converged=true;
            runMore=false;
        end
        if(C.iter>=C.maxIter) % check for max iterations
            runMore=false;
            C.exitStatus=[C.exitStatus 'Maximum number of iterations reached. '];
        end
        if(displayProgress) % display convergence progress
            if(mod(C.iter,10)==1)
                displayHeader();
            end
            fm=displayConvergence();
            if(fm)
                disp('---------- lower bound decrease?!!! ----------')
            end
        end
        if(C.iter<2 && C.iter<C.maxIter) % run at least 2 iterations unless specifi
            runMore=true;
            C.exitStatus='';
        end
    end
    %% estimates (only last iteration)
    if(~runMore)        
        % global light-weight estimates (always)
        wa0=sum(W.M.wa,2);
        W.est.aMean=W.M.wa(:,1)./wa0;
        W.est.aMode=(W.M.wa(:,1)-1)./(wa0-2);
        W.est.aVar=W.M.wa(:,1).*W.M.wa(:,2)./(wa0.^2.*(1+wa0));
        
        wB0=sum(W.M.wB,2)*ones(1,W.N);
        eyeB=1-eye(W.N);
        W.est.Bmean=W.M.wB./wB0;
        W.est.Bmode=(W.M.wB-1+eye(W.N))./(wB0-W.N+1);
        W.est.Bvar=W.M.wB.*(wB0.*eyeB-W.M.wB)./(wB0.^2.*(1+wB0));
        B2=W.M.wB.*(eyeB+W.M.wB)./(wB0.*(1+wB0));  % <Bjk^2>
        
        W.est.Amean=diag(W.M.wa(:,2)./sum(W.M.wa,2))...
            +(W.M.wa(:,1)./sum(W.M.wa,2)./sum(W.M.wB,2))*ones(1,W.N).*W.M.wB;
        %W.est.Amode : have not figured that one out yet (ML 2014-05.02)
        W.est.Astd=diag(W.est.aVar)...
            +(W.est.aVar*ones(1,W.N)).*B2...
            +(W.est.aMean.^2*ones(1,W.N)).*W.est.Bvar;
        
        W.est.dwellMean=1./W.est.aMean;
        W.est.dwellMode=wa0./(1+W.M.wa(:,1));
        clear wB0 eyeB B2 wa0

        % emission parameters %%% remains to do
        [~,~,~,~,~,lambda,W.est.sig2Mean]=...
            VB5_ag_vl_expectations_beta0(W);
        W.est.DMean=lambda/2/W.param.timestep;
        %%%W.est.DMode=W.M.cg/2./(W.M.ng+1)/W.param.timestep;
        
        % occupation
        W.est.Ttot=sum(W.Es.pst,1);
        W.est.Ptot=W.est.Ttot/sum(W.est.Ttot);

        W.Fterms.Fterms=[ W.Fterms.lnZQ+W.Fterms.lnZq+W.Fterms.lnZz -sum(KL_a) -sum(KL_B) -sum(KL_pi) -sum(KL_lv)];
        W.Fterms.FtermsNames='[lnZQ+lnZq+lnZz -sum(KL_a) -sum(KL_B) -sum(KL_pi) -sum(KL_lv)]';

        W.est.Ps=W.est.Ttot/sum(W.est.Ttot);
        
        %% potentially demanding estimates (only if asked)
        if(do_estimates)
            % extract trajectory estimates
            Wviterbi=uint8(HMM_multiViterbi_log(lnQ,lnH,trjEnd)); % Viterbi path
            [~,WsMaxP]=max(W.Es.pst,[],2);
            
            for kk=1:length(trjStart)
                W.est2.pst{kk}    =      W.Es.pst(trjStart(kk):trjEnd(kk),:);
                W.est2.H{kk}      =             H(trjStart(kk):trjEnd(kk),:);
                W.est2.lnH{kk}    =           lnH(trjStart(kk):trjEnd(kk),:);
                W.est2.lnHMax{kk} =        lnHMax(trjStart(kk):trjEnd(kk),:);
                W.est2.sMaxP{kk}  =uint8(  WsMaxP(trjStart(kk):trjEnd(kk)));
                W.est2.viterbi{kk}=uint8(Wviterbi(trjStart(kk):trjEnd(kk)));
                %W.est2.start=trjstarts;
                %W.est2.end=trjEnds;
            end
            clear Wviterbi WsMaxP
            
            W.est2.Ts=zeros(length(trjStart),N);
            W.est2.Ps=zeros(length(trjStart),N);
            for m=1:length(trjStart)
                W.est2.Ts(m,:)=sum(W.Es.pst(trjStart(m):trjEnd(m),:),1); % time spent in each state
                W.est2.Ps(m,:)=W.est2.Ts(m,:)/sum(W.est2.Ts(m,:));
            end
            
            try
                W.est.lnAmean=logm(W.est.Amean);
                W.est.lnA_error='none';
            catch me
                W.est.lnAmean=0*W.est.Amean;
                W.est.lnA_error=me;
            end
        end
    end
end
%% exit message
if(displayExit) % display exit message
    displayHeader();
    displayConvergence();
    fprintf('%s \n',C.exitStatus)
end
%% slim down the model, on request, by deleting the bulky E-field
%%% ML: maybe some other bulky fields should go as well?
if(do_slim)
    W=rmfield(W,{'E','Etrj'});
    W.est=rmfield(W.est,{'Ts','Ps'});
end
%% auxiliary functions
    function dFminus=displayConvergence()
        fprintf('%02d % 5d % 0.2e ',[W.N C.iter C.dFrel]);
        fprintf('%0.2e  %s ',C.dPrel,C.limitPar)
        %fprintf(', dFterms :')
        %fprintf('% 0.2e ',[W.Fterms.Fterms-Wm1.Fterms.Fterms]);
        fprintf(', F : %0.3e \n',W.F)
        dFminus=(C.dFrel<-1e-11);
        if(dFminus)
            %keyboard
        end
    end
    function displayHeader()
        disp(' N     it   dF/|F|   d<p>/p   p ')
    end
end

