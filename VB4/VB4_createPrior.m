function W=VB4_createPrior(runinput,N)
% W=V3_createPrior(runinput,N)
%
% Creates a model structure W with N states, and prior distributions
% according to runinput, either a runinputfile or an options struct. By
% default no aggregated states are generated. This file calls
% VB3_createPrior for creating those parts of the model that are common to
% both types of analysis. 
% 
% Note: VB3 and VB4 currently uses different conventions for parameterizing
% the diffusion constant.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB4_createPrior.m, model initialization in the vbSPT package
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
%% Parse input and generate VB3 prior
% if an existing file, generate options structure
if(ischar(runinput) && exist(runinput,'file'))
    runinputfile = runinput;
    opt=VB3_getOptions(runinputfile);
    % if an option struct, read in the runinputfilename
elseif(isstruct(runinput))
    opt=runinput;
else
    error(['Not a valid input, aborting VB4_createPrior']);
end
W=VB3_createPrior(opt,N);
W.PM=rmfield(W.PM,{'n','c'}); % obsolete fields from VB3 model

%% start of actual code
W.param.blur_R=opt.blur_R;
W.param.blur_tau=opt.blur_tau;
W.param.blur_beta=opt.blur_tau*(1-opt.blur_tau)-opt.blur_R;
W.param.timestep=opt.timestep;

%% default prior type choices for VB4
if(~isfield(opt,'prior_type_locErr'))
    warning('VB4_createPrior: prior_type_locErr not specified. Using default : mean_strengt')
    prior_type_D='mean_strength';
end
%% localization error prior
if(strcmp(opt.prior_type_locErr,'mean_strength'))
    [W.PM.na,W.PM.ca]=priorLocErr_mean_strength(opt,N);
else
    error(['VB4_createPrior: did not recognize prior_type_locErr : ' opt.prior_type_locErr])
end
%% diffusion constant prior
if(strcmp(opt.prior_type_D,'mean_strength'))
    [W.PM.ng,W.PM.cg]=priorD_mean_strength(opt,N);
else
    error(['VB4_createPrior: did not recognize prior_type_D : ' opt.prior_type_D])
end
end
%% slightly complicated prior choices

function [m,h]=priorLocErr_mean_strength(opt,N)
% Default dynamic localization error prior, specified using the mean value
% and total strength per state.

E0=opt.prior_locErr;                     % prior mean of dynamic localization error
En=opt.prior_locErrStrength;             % prior strength of dynamic localization error

% each emission variable gets same strength independent of model size
m=En*ones(1,N)+1;
h=E0^2*(m-1);
end

function [n,c]=priorD_mean_strength(opt,N)
% Default diffusion constant prior, specified using the mean value and
% total strength per state. VB4 has a different definition of gamma than
% VB3, to get rid of several confusing factors of two.

timestep=opt.timestep;              % sampling time step
D0=opt.prior_D;                     % prior diffusion constant
Dn=opt.prior_Dstrength;             % strength of diffusion constant prior

% each emission variable gets same strength independent of model size
n=Dn*ones(1,N)+1;
c=2*D0*timestep*(n-1); % match <gamma>=<1/2Ddt> 

end
