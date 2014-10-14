function W=VB5_createPrior(runinput,N)
% W=VB5_createPrior(runinput,N)
%
% Creates a model structure W with N states, and prior distributions
% according to runinput, either a runinputfile or an options struct. By
% default no aggregated states are generated. This file calls
% VB3_createPrior for creating those parts of the model that are common to
% both types of analysis.
%
% Note: VB3 and VB5 currently uses different conventions for parameterizing
% the diffusion constant.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB5_createPrior.m, model initialization in the vbSPT package
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
    error('Not a valid input, aborting VB5_createPrior.');
end
W=VB3_createPrior(opt,N);
W.PM=rmfield(W.PM,{'n','c'}); % obsolete fields from VB3 model

%% start of actual code
W.param.blur_R=opt.blur_R;
W.param.blur_tau=opt.blur_tau;
%W.param.blur_beta=opt.blur_tau*(1-opt.blur_tau)-opt.blur_R;
W.param.timestep=opt.timestep;

%% localization error prior
if(~isfield(opt,'prior_type_locErr'))
    error('VB5_createPrior: prior_type_locErr not specified.')
elseif(strcmp(opt.prior_type_locErr,'rms_strength'))
    sigma2=opt.prior_locErr^2;
    strength=opt.prior_locErrStrength;
    [W.PM.nv,W.PM.cv]=inverseGamma_mean_strength(sigma2,strength,N);
elseif(strcmp(opt.prior_type_locErr,'mean_strength'))
    sigma=opt.prior_locErr;
    strength=opt.prior_locErrStrength;
    [W.PM.nv,W.PM.cv]=inverseGamma_meanSqrt_strength(sigma,strength,N);
else
    error(['VB5_createPrior: did not recognize prior_type_locErr : ' opt.prior_type_locErr])
end
%% diffusion constant prior
if(~isfield(opt,'prior_type_D'))
    error('VB5_createPrior: prior_type_D not specified.')
elseif(strcmp(opt.prior_type_D,'mean_strength'))
    lambda=2*opt.prior_D*W.param.timestep;
    strength=opt.prior_Dstrength;
    [W.PM.nl,W.PM.cl]=inverseGamma_mean_strength(lambda,strength,N);
else
    error(['VB5_createPrior: did not recognize prior_type_D : ' opt.prior_type_D])
end
end
%% slightly complicated prior choices
function [n,c]=inverseGamma_mean_strength(mean,strength,N)
% Default dynamic localization error prior, specified using the mean value
% and total strength per state.
if(strength<=1)
    error('inverseGamma_mean_strength: prior strength must be > 1 for mean strength option.')
end
% each emission variable gets same strength independent of model size
n=strength*ones(1,N);
c=mean*(strength-1)*ones(1,N);
end

function [n,c]=inverseGamma_meanSqrt_strength(mean,strength,N)
% Default dynamic localization error prior, specified using the mean value
% and total strength per state.
if(strength<=-0.5)
    error('inverseGamma_meansqrt_strength: prior strength must be > -0.5 for meansqrt strength option.')
end
% each emission variable gets same strength independent of model size
n=strength*ones(1,N);
if(strength>1000)
    warning('inverseGamma_meansqrt_strength can be numerically unstable for strength > 1000')
end
c=mean^2*exp(2*(gammaln(strength)-gammaln(strength-0.5)))*ones(1,N);
end
