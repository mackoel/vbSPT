function dat=VB4_preprocess(runinput,dim)
% dat=VB4_preprocess(runinput,dim)
%
% assemble single particle diffusion data in a form for fast EM iterations
%
% runinput  : a VB3/VBSPT runinput file, or a VB3/VBSPT runinput structure,
%             or a cell vector of diffusion trajectories
% optional parameters with cell vector input (ignored otherwise):
% dim       : optional specification of data dimension when input is cell
%             vector; only the first dim columns of each trajectory will be
%             analysed. Default: use all columns. If options structure or
%             runiput file is given, this argument is gnored.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB4_preprocess, data preprocessor in the vbSPT package
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
%% parse input
% if an existing file, generate options structure
if(ischar(runinput) && exist(runinput, 'file')==2)
    runinputfile = runinput;
    opt=VB3_getOptions(runinputfile);
    disp(['Read runinput file ' runinputfile])
    X=VB3_readData(opt);
    dim=opt.dim;
elseif(isstruct(runinput)) % if already an options structure
    opt=runinput;
    %runinputfile=opt.runinputfile;
    X=VB3_readData(opt);
    dim=opt.dim;
elseif(iscell(runinput))
    X=runinput;
    if(~exist('dim','var'))
       dim=size(X{1},2);
    end
else
    error(['Not a valid input, aborting SPT_preprocess']);
end
clear opt runinput runinputfile
%% assemble output structure
dat=struct;
dat.dim=dim;

% count output size
T=zeros(size(X));
for k=1:length(X)
    T(k)=size(X{k},1);
end
if(~isempty(find(T<2,1)))
   warning('VB4_preprocess: data contains traces with no steps.')
end
% data stacking: pack a zero-row between every trajectory to match sizes of
% data x(t) and diffusive path y(t).
dat.T=T;
dat.x=zeros(sum(T),dim);
dat.one=zeros(1,length(X),'double');
dat.end =zeros(1,length(X),'double');

ind=1;
for k=1:length(X)
    x=X{k}(:,1:dim);
    Tx=size(x,1);
    dat.one(k)=ind;    
    dat.end(k)  =ind+Tx-1;
    ind=ind+Tx;
    dat.x(dat.one(k):dat.end(k),1:dim)=x;
end




