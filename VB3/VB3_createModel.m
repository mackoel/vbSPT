function M=VB3_createModel(D,A,p0,W0,dt,strength)
% function M=VB3_createModel(D,A,p0,W0,dt,strength)
% 
% create a VB3 model from parameters and strength. The function either
% returns an M-field, or (if an input model W0 is given), adds the M-field
% to an input model. The model parameters are chosen so that the
% vatiational mean values coincide with the input parameters. This is good
% when constructing initial guesses for synthetic data with known
% parameters.
%
% D  : vector of diffusion constants * timestep
% A  : transition matrix, needs to be properly normalized etc
% p0 : initial state probability (default = approx. stationary state of A)
% strength  : strength parameter, corresponding to the number of counts that go
%      into each parameter distribution (default = 1e4). 
% W0 : optional starting model. Only the PM and dim fields of this
%      model is inherited, if present. If not given, only the M-field
%      corresponding to the other input parameters are returned.
% dt : data timestep (default=1, which means that one can omit dt and input
%      D*dt for the diffusion constants directly)

% M.L. 2012-07-05

%% check input parameters and size compatibility
if(exist('dt','var') && ~isempty(dt) && length(dt)==1 && dt>0)
   Ddt=D*dt;
else
    if(length(dt)~=1 || dt<=0)
        error('VB3_createModel: dt must be a positive scalar')
    end
   Ddt=D;
end
%clear dt D;

lambda=sort(eig(A));
if(abs(lambda(end)-1)>1e-15) % then there is no proper stationary state
    error('VB3_createModel: A does not seem to be a proper transition matrix')
end

if(~exist('p0','var') || ~isempty(p0)) % then use stationary state of A
    % compute an equilibration time from second largest eigenvalue
    % at this point, one knows that lambda(end-1)<1
    if(abs(real(lambda(end-1)))<1)
        neq=-1/log(abs(real(lambda(end-1))));
        p0=A^ceil(100*neq);
        p0=p0(1,:);
    else
        error('VB3_createModel: cannot compute equilibration time for transition matrix')
    end
end
%clear lambda

N_Ddt=length(Ddt);
N_A=size(A);
N_p0=length(p0);
if( N_Ddt==N_A(1) && N_A(1)==N_A(2) && N_A(2)==N_p0)
else
    error('VB3_createModel: incompatible input parameter sizes!')
end
N=N_Ddt;
clear N_Ddt N_A N_p0

if(~exist('strength','var') || isempty(strength))
    strength=1e4;
end
if(length(strength)~=1 || strength<=1)
    error('VB3_createModel: strength must be a scalar > 1')
end

%% compute M field
M.n=ones(1,N)*strength;
M.c=4*Ddt*(strength-1);
M.wA=A*strength;
M.wPi=p0*strength;
%clear strength Ddt D A p0

% assemble result structure
if(exist('W0','var')) % thgen we shall return an output model
    Wout=struct;
    Wout.M=M;
    if(isfield(W0,'PM'))
        Wout.PM=W0.PM;
    end
    if(isfield(W0,'dim'))
        Wout.dim=W0.dim;
    end
    if(isfield(W0,'N') && W0.N~=N)
        error('VB3_createModel: input model size incompatible with parameters')
    end    
    M=Wout;
end