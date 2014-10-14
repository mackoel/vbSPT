function M=VB5_createModel(Ddt,locErr,A,p0,strength,W0)
% function M=VB5_createModel(Ddt,locErr,A,p0,strength,W0)
% 
% Adds a new M-field M to a VB5 model W0, based on given parameters and
% strength. Model variational parameters are chosen so that the vatiational
% mean values coincide with the input parameters.
%
% Ddt: vector of diffusion constants * timestep
% locErr: vector of localization errors for each state
% A  : transition matrix, needs to be properly normalized etc
% p0 : initial state probability (default = approx. stationary state of A)
% strength : strength parameter, corresponding to the number of counts
%            that go into each parameter distribution (default = 1e4). 
% W0 : Starting model, for example created by VB5_createPrior, with an
%      existing PM field of matching dimensionality.
%
% M.L. 2014-09-25

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB5_createModel, creates a VB5 model from parameters and strengths.
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
    error('VB5_createModdt,el: A does not seem to be a proper transition matrix')
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
            error('VB5_createModel: cannot compute equilibration time for transition matrix')
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
    error('VB5_createModel: incompatible input parameter sizes!')
end
clear N_Ddt N_A N_p0

if(~exist('strength','var') || isempty(strength))
    strength=1e4;
end
if(length(strength)~=1 || strength<=1)
    error('VB5_createModel: strength must be a scalar > 1')
end

%% compute M field
M.nl=ones(1,N)*strength;
M.cl=Ddt*(strength-1);

% if this method is to work without an input model, then a beta value must
% be supplied instead.
tau=W0.param.blur_tau;
beta=tau*(1-tau)-W0.param.blur_R;
M.na=ones(1,N)*strength;
M.ca=M.na.*(2*Ddt*beta+locErr.^2);

wA=A*strength;
M.wa=[sum(wA,2)-diag(wA) diag(wA)];
M.wB=wA-diag(diag(wA));

M.wPi=p0*strength;

% check size match, and add M-field to W0 if OK. These checks are somewhat
% redundant, but will be useful if the method is extended to do without an
% input model.
if(~exist('W0','var'))
    return
elseif(~isstruct(W0) || ~isfield(W0,'PM'))
    error('VB5_createModel: W0 is not a struct, or PM field is missing.')
else
    Mname=fieldnames(M);
    Pname=Mname;
    for k=1:length(Mname)
        if(strcmp(Mname{k},'na'))
            Pname{k}='nv';
        elseif(strcmp(Mname{k},'ca'))
            Pname{k}='cv';
        end
            
        if(~isfield(W0.PM,Pname{k}))
            error(['VB5_createModel: field ' Pname{k} ' missing from W0.PM'])
        else
            [a0,b0]=size(W0.PM.(Pname{k}));
            [a1,b1]=size(M.(Mname{k}));
            if(a0~=a1 || b0~=b1)
                error(['VB5_createModel: dimensionality of field W0.PM.' Pname{k} ' does not match input parameters'])
            end
        end
    end
    M.SA=1:N;
    W0.M=M;
    W0.N=N;
    M=W0;
end
