function [x,s,y]=VB5_diffusiveHMM_blur_detach(p0,A,pE,Ddt,R,tau,sigErr,dim,T,numTrj)
% [x,s,y]=VB5_diffusiveHMM_blur_detach(p0,A,pE,Ddt,R,tau,sigErr,dim,T,numTrj)
%
% simulate hidden states s(t), true positions y(t), and measured positions 
% x(t) for a diffusion that switches between different diffusion constant 
% according to a Markov process. Motional blur and localization errors are 
% added according to the Michalet&Berglund model.
% Cell vectors are returned if more than one trajectory is simulated.
%
% p0 : initial state distribution. Default: stationary state of A. p0 is
%      automatically normalized.
% A  : transition matrix, using the convention 
%      A(i,j)=p(s_t=j|s_{t-1}=i), so the evolution of the hidden states is
%      given by p(t)=p(t-1)*A, where p(t) is a row vector
% pE : trajectory end probability vector, pE(j)=p(t is last point|s_t=j).
%      Default 0, and a scalar is interpreted as same for all states.
%      Minimum trajectory length is 2.
%
% R  : motional blur factor 0<R<1/4
% tau: mean value of illumination density, in units of the timestep dt
% sigErr: standard deviation of state-dependent Gaussian localization error
%
% dim: spatial dimension of the output data (default 2)
%
% T  : maximum trajectory length(s). Only the first M entries are used, and
% if numTrj>length(T), then the entries are cycled. Default: 100.

% numTrj : number of trajectories to simulate. Default: 1.
% Ddt: diffusion*timestep, Ddt(j)=D(j)*dt
%
% ML 2014-08-21, bmelinden@gmail.com

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB5_diffusiveHMM_blur_detach.m, simulate d-dimensional diffusion with 
% multiple diffusion constants, Michalet-Berglund blur model, and state-
% dependent detachments, part of the vbSPT package
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

%% parameter handling
[N,NN]=size(A);
if(N~=NN)
    error('VB5_diffusiveHMM_blur_detach requires square transition matrix')
end
clear NN

if( isempty(p0) || length(p0)==1 )
    A0=A^10000;
    p0=A0(1,:);
end
p0=p0/sum(p0);

if(isempty(pE))
    pE=zeros(1,N);
elseif(length(pE)==1)
    pE=pE*ones(1,N);
end
if(R<0 || R>0.25)
    error('R-values must be in the range 0<R<0.25 (R<1/6 for constant illumination).')
end

if(tau<0 || tau>1)
    error('Illumination mean time tau must be 0<tau<1.')
end
if(R > tau*(1-tau))
   error('R > tau*(1-tau), which is unphysical.')
end
if(length(sigErr)==1)
    sigErr=sigErr*ones(1,N);
end
if(~exist('dim','var') || isempty(dim)); dim=2; end
if(~exist('T','var') || isempty(T)); T=100; end
if(~exist('numTrj','var') || isempty(numTrj)); numTrj=1; end

Ddt=reshape(Ddt,length(Ddt),1);
sigErr=reshape(sigErr,length(sigErr),1);


%% initialization
x=cell(1,numTrj);
s=cell(1,numTrj);
y=cell(1,numTrj);
%% start simulation
NT=length(T);
cumA=cumsum(A,2);
cum0=cumsum(p0);
for m=1:numTrj
    % initialize trajectory
    Ttrj=T(1+mod(m-1,NT)); % max length of this trajectory
    if(Ttrj<inf)
        S=zeros(Ttrj,1);
    else
        S=zeros(10,1);
    end
    % hidden state trajectory
    S(1)=find(rand<cum0,1); % initial state
    for t=2:Ttrj
        ras=rand;
        S(t)=find(ras<cumA(S(t-1),:),1);
        if(rand<pE(S(t))) % then terminate trajectory here
            S=S(1:t);
            break
        end
    end
    Ttrj=length(S);
    
    % pure diffusion-trace
    dY=zeros(Ttrj,dim);
    for k=1:dim
        dY(:,k)=randn(Ttrj,1).*sqrt(2*Ddt(S));
    end
    Y=[zeros(1,dim); cumsum(dY,1)];
    % add motional blur and localization noise
    X=Y(1:Ttrj,:)*(1-tau)+Y(2:end,:)*tau;
    for k=1:dim
        X(:,k)=X(:,k)+randn(Ttrj,1).*...
            sqrt(2*Ddt(S)*(tau*(1-tau)-R)+sigErr(S).^2);
    end
    if(numTrj==1)
        x=X;s=S;y=Y;
    else
        x{m}=X;s{m}=S;y{m}=Y;
    end
end
