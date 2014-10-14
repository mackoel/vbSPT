function w=VB5_removeState(w,s,opt)
% w=VB5_removeState(w,s,opt)
% create a new model parameter field M by removing state s from model w.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VB5_removeState.m, removes state from models in the vbSPT package
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
%  If you modify this Program, or any covered work, by linking or combining
%  it with Matlab or any Matlab toolbox, the licensors of this Program
%  grant you additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code
M=struct;

sk=[1:s-1 s+1:w.N]; % states to keep
M.wPi=w.M.wPi(sk);
M.nl  =  w.M.nl(sk);
M.cl  =  w.M.cl(sk);
M.na  =  w.M.na(sk);
M.ca  =  w.M.ca(sk);
M.wa =  w.M.wa(sk,:);
M.SA =  w.M.SA(sk);
M.SA(M.SA>s)=M.SA(M.SA>s)-1; % decrease aggregation indices

% transfer observed transitions
wB =w.M.wB -w.PM.wB;
% try to compensate for transitions that went via the removed state
toS=wB(sk,s);
frS=wB(s,sk);
M.wB=wB(sk,sk)+w.PM.wB(sk,sk)...
    +(toS*frS/sum(frS)+toS*frS/sum(toS)).*(1-eye(length(sk)));

if(size(M.wa,1)==1) 
    % B and a are not defined for only one state, but M.wa is still kept for
    % numerical ease (and is taken care of elswhere).
    M.wB=0;
end

Etrj=w.Etrj;
w=VB5_createPrior(opt,w.N-1);
w.M=M;
w.Etrj=Etrj;
