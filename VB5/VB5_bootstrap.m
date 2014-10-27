function [wbs,Wmean,Wstd]=VB5_bootstrap(W,dat,opt,NB,varargin)
% [wbs,Wmean,Wstd]=VB5_bootstrap(W,dat,opt,NB,wbs0)
%
% Run NB bootstrap fits starting from model W and data X. Each bootstrap
% fit resamples the data with replacement, and then converges the VB5 model
% using W as a starting point. No other model sizes are attempted, so if
% the data set is small, or if some states are only present in very few
% trajectories, strange things may happen.
%
% No path estimates are done, even if opt.pathestimates=true
% This function uses a parfor loop, but does not open or close the
% matlabpool by itself. 
%
% W     : Initial guess for bootstrap fit. For good results, use the best
%         model for the inriginal data set X. 
% X     : diffusion data 
% opt   : VB5 parameter object ( obj=VB5_getOptions(runinputfile) ).
% NB    : optional number of bootstrap iterations to run. Default: 100.
% wbs0  : wbs struct array, containing bootstrap indices in wBS.ind. This
%         is used to duplicate a bootstrap analysis on a different starting
%         model.
%
% wbs   : vector of bootstrapped models. To save some space, only three fields are stored. 
%         wbs(k).M, wbm(k).est (from the bootstrapped models), and
%         wbs(k).ind, the resampling indices which can be used to recreate
%         the resampled data, as Y={X{wbs(k).ind}};, and from there the
%         rest of the model.
% Wmean, Wstd: bootstrap average and standard deviation models, taken
%              elements-wise for some chosen fields. Only M and est fields
%              are analyzed.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB5_bootstrap.m, low-level bootstrap analysis in the VB5 suite of vbSPT
% =========================================================================
% 
% Copyright (C) 2014 Martin Lindén, bmelinden@gmail.com
%
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
%% parameters
hasindices=false;
if(nargin>4)
    wbs0=varargin{1};
    NB=length(wbs0);
    hasindices=true;
else
    wbs0=struct('ind',cell(1,NB)); % trick to avoid complaints from parfor
end
if(~exist('NB','var') || isempty(NB)); NB=100;end
%% bootstrap fits
L=length(dat.one);
tic
%%%parfor k=1:NB
for k=1:NB
    %tic
    if(hasindices)
        ind=wbs0(k).ind;
    else
        ind=sort(ceil(L*rand(1,L))); % sample with replacement
    end
    %disp(['bootstrap iter ' int2str(k) ' 1'])
    % produce resampled data set
    datY=dat;
    ww=W;
    
    % resample trajectories and input guess 
    tic
    EtrjOne=1;
    EtrjEnd=0;
    datOne=1;
    datEnd=0;
    for tt=1:length(datY.one)
        TT=dat.T(ind(tt));
        
        % resample old data
        datY.T(tt)=TT;
        datY.one(tt)=datOne;
        datEnd=datOne+TT-1;
        datY.end(tt)=datEnd;
        datOne=datEnd+1;
        x0= dat.one(ind(tt)):dat.end(ind(tt));
        x1=datY.one(tt):datY.end(tt);
        datY.x(x1,:)=dat.x(x0,:);
        
        % resample old input guess
        ww.Etrj.one(tt)=EtrjOne;
        EtrjEnd=EtrjOne+TT;
        ww.Etrj.end(tt)=EtrjEnd;
        EtrjOne=EtrjEnd+1;
        
        mu0= W.Etrj.one(ind(tt)):W.Etrj.end(ind(tt));
        mu1=ww.Etrj.one(tt):ww.Etrj.end(tt);
        
        ww.Etrj.mu(mu1,:)      =W.Etrj.mu(mu0,:);
        ww.Etrj.CovDiag0(mu1)=W.Etrj.CovDiag0(mu0);
        ww.Etrj.CovDiag1(mu1(1:end-1))=W.Etrj.CovDiag1(mu0(1:end-1));
    end
    datY.x=datY.x(1:datEnd,:);    
    ww.Etrj.mu=ww.Etrj.mu(1:EtrjEnd,:);
    ww.Etrj.CovDiag0=ww.Etrj.CovDiag0(1:EtrjEnd);
    ww.Etrj.CovDiag1=ww.Etrj.CovDiag1(1:EtrjEnd-1);
    disp(['Bootstrap resampled data and model : ' num2str(toc) ' s.'])
    ww=rmfield(ww,{'Es','Epar','E','est'});
    
    %disp(['bootstrap iter ' int2str(k) ' 2'])
    save BSsnapshot.mat %% for debugging
    ww=opt.VBEMfunction(ww,datY,'outputLevel',2,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar,'slim');
    %disp(['bootstrap iter ' int2str(k) ' 3'])    
    wbs(k).M=ww.M;
    wbs(k).est=ww.est;
    wbs(k).F=ww.F;    
    wbs(k).ind=ind;
    disp(['bootstrap iter ' int2str(k) ' finished'])
end
disp(['bootstrap : ' int2str(NB) ' resampling in ' num2str(toc) ' s.'])
%% default bootstrap analysis
Wmean=struct;
Wstd=struct;
f={'M' 'est'}; % fields on which to compute bootstrap statistics

for m=1:length(f)
   g=fieldnames(wbs(1).(f{m}));
   for j=1:length(g)
       [a,b]=size(wbs(1).(f{m}).(g{j}));
       bsnum=zeros(NB,a*b);
       
       % construct table
       for k=1:NB
            bsnum(k,:)=wbs(k).(f{m}).(g{j})(1:end);
       end
       % compute mean and std
       Wmean.(f{m}).(g{j})       =zeros(a,b);
       Wmean.(f{m}).(g{j})(1:end)=mean(bsnum,1);       
       Wstd.(f{m}).(g{j})       =zeros(a,b);
       Wstd.(f{m}).(g{j})(1:end)=std(bsnum,[],1);
   end
end
Wmean.F=mean([wbs(1:end).F]);
Wstd.F=std([wbs(1:end).F]);


