function M=VB4_createModel(Ddt,locErr,A,p0,strength,W0)
% function M=VB4_createModel(Ddt,locErr,A,p0,strength,W0)
% 
% Constructs a new M-field M to a VB4 model, based on given parameters and
% strength. Model variational parameters are chosen so that the vatiational
% mean values coincide with the input parameters. If an input model W0 is
% given, the M field is instead added to that model after some
% compatibility checks, and the whole model struct returned.
%
% Ddt: vector of diffusion constants * timestep
% locErr: vector of localization errors for each state
% A  : transition matrix, needs to be properly normalized etc
% p0 : initial state probability (default = approx. stationary state of A)
% strength : strength parameter, corresponding to the number of counts
%            that go into each parameter distribution (default = 1e4). 
% W0 : Optional starting model, for example created by VB4_createPrior,
%      with an existing PM field of matching dimensionality.
%


%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB4_createModel, creates a VB4 model from parameters and strengths.
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

%% check M-field input parameters and size compatibility    
% check transition matrix
lambda=sort(eig(A));
if(abs(lambda(end)-1)>1e-15) % then there is no proper stationary state
    error('VB4_createModdt,el: A does not seem to be a proper transition matrix')
end
if(~exist('p0','var') || ~isempty(p0)) % then use stationary state of A
    % compute an equilibration time from second largest eigenvalue
    % at this point, one knows that lambda(end-1)<1
    if(size(A,1)>1)
        if(abs(real(lambda(end-1)))<1)
            neq=-1/log(abs(real(lambda(end-1))));
            p0=A^ceil(100*neq);
            p0=p0(1,:);
        else
            error('VB4_createModel: cannot compute equilibration time for transition matrix')
        end
    else
        p0=1;
    end
end
clear lambda

N_Ddt=length(Ddt);
N_A=size(A);
N_p0=length(p0);
if( N_Ddt==N_A(1) && N_A(1)==N_A(2) && N_A(2)==N_p0)
    N=N_Ddt;
else
    error('VB4_createModel: incompatible input parameter sizes!')
end
clear N_Ddt N_A N_p0

if(~exist('strength','var') || isempty(strength))
    strength=1e4;
end
if(length(strength)~=1 || strength<=1)
    error('VB4_createModel: strength must be a scalar > 1')
end

%% compute M field
M.ng=ones(1,N)*strength;
M.cg=Ddt*(strength-1);

M.na=ones(1,N)*strength;
alpha=1./locErr.^2;
M.ca=M.na./alpha;

wA=A*strength;
M.wa=[sum(wA,2)-diag(wA) diag(wA)];
M.wB=wA-diag(diag(wA));

M.wPi=p0*strength;

% check size match, and add M-field to W0 if OK
if(~exist('W0','var'))
    return
elseif(~isstruct(W0) || ~isfield(W0,'PM'))
    error('VB4_createModel: W0 is not a struct, or PM field is missing.')
else
    fn=fieldnames(M);
    for k=1:length(fn)
        if(~isfield(W0.PM,fn{k}))
            error(['VB4_createModel: field ' fn{k} ' missing from W0.PM'])
        else
            [a0,b0]=size(W0.PM.(fn{k}));
            [a1,b1]=size(M.(fn{k}));
            if(a0~=a1 || b0~=b1)
                error(['VB4_createModel: dimensionality of field ' fn{k} ' does not match input parameters'])
            end
        end
    end
    M.SA=1:N;
    W0.M=M;
    W0.N=N;
    M=W0;
end