function W0=VB5_trjInit_dat(dat,W0,muStd,muCov)
% W0=VB5_trjInit_dat(dat,W0,muStd,muCov)
% 
% Creates a trajectory field Etrj based on the VB5 data dat, that can be
% used as an initial guess. If an input model is given, a model is returned
% with the new Etrj field added. The data object dat can be created with
% VB4_preprocess.
%
% muStd : initial standard deviation of hidden trajectory (default 0).
% muCov : initial nearest-neightbour covariance of hidden trajectory,
%         <mu(t)mu(t+1)> (default 0).

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB5_trjInit_dat, creates an initial VB4 hidden trajectory. Part of the
% vbSPT suite.
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
%% Input parameters for optional Etrj field
if(~exist('muStd','var') || isempty(muStd))
    muStd=0;
end
if(~exist('muCov','var') || isempty(muCov))
    muCov=0;
end
%% create Etrj field
Etrj=struct;
dim=dat.dim; 
% indices to hidden trajectory start- and end-points
Etrj.one=dat.one+(0:length(dat.one)-1);
Etrj.end=dat.end+(1:length(dat.one)  );

% create mean values
Etrj.mu=zeros(sum(dat.T+1),dim);
for k=1:length(Etrj.one)
   muk=dat.x(dat.one(k):dat.end(k),:);
   for d=1:dim
      muk(:,d)=smooth(muk(:,d),3); 
   end
   Etrj.mu(Etrj.one(k):Etrj.end(k)-1,:)=muk;
end
Etrj.mu(Etrj.end,:)=dat.x(dat.end,:);

Etrj.CovDiag0=muStd^2*ones(size(Etrj.mu,1),1);
Etrj.CovDiag1=muCov*ones(size(Etrj.CovDiag0));
Etrj.CovDiag1(Etrj.end)=0;

%% add the Etrj field to W0 (or return Etrj of W0 does not exist).
if(~exist('W0','var')|| isempty(W0))    
    W0=Etrj;
else
    if(isfield(W0,'dim') && W0.dim~=dim)
        warning('VB5_createTrajectory: model dimensionality changed.')
    end    
    W0.Etrj=Etrj;
    W0.dim=dim;
end
