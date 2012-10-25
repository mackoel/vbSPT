function est=VB3_parameterEstimates(M)
% est=VB3_parameterEstimates(M)
% generate some estimates based on a M-field. Good for interpreting various
% 'true' models, e.g., those generated by knowing the hidden states etc.
% M : input M field of some model
% est : those parameter estimates that can be done based on the M field
% alone, same as done by VB3_VBEMiterator

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_parameterEstimates, generate some estimates based on a M-field
% =========================================================================
% 
% Copyright (C) 2012 Martin Lind??n and Fredrik Persson
% 
% E-mail: bmelinden@gmail.com, freddie.persson@gmail.com
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



N=length(M.wPi);
% global light-weight estimates (always)
est.Amean=zeros(size(M.wA));
est.Astd=zeros(size(M.wA));
est.Amode=zeros(size(M.wA));
a0=sum(M.wA,2);
for k=1:N
    est.Amean(k,:)=M.wA(k,:)/a0(k);
    est.Astd(k,:)=sqrt((M.wA(k,:).*(a0(k)-M.wA(k,:)))./(a0(k)^2*(1+a0(k))));
    est.Amode(k,:)=(M.wA(k,:)-1)/(a0(k)-N);
end

est.dwellMean=1./(1-diag(est.Amean));
est.dwellMode=1./(1-diag(est.Amode));
% emission parameters
est.gMean=M.n./M.c;
est.gMode=(M.n-1)./M.c;
est.gStd=sqrt(M.n./M.c.^2); % sqrt(Var(g))
est.DdtMean=M.c/4./(M.n-1);
est.DdtMode=M.c/4./(M.n+1);
est.Ddtstd=M.c/4./(M.n-1)./sqrt(M.n-2);
