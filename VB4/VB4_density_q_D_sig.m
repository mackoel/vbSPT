% function to plot the variational posterior density q(D,\sigma)
% [D,sig,rho,ln_rho]=VB4_density_q_D_sig(W,D_in,sig_in,DSbins)

function [D,sig,rho,ln_rho]=VB4_density_q_D_sig(W,D_in,sig_in,DSbins)

% parse input
if(exist('DSbins','var') && length(D_in)==2 && length(sig_in)==2 )
   D_in=linspace(D_in(1),D_in(2),DSbins(1));
   sig_in=linspace(sig_in(1),sig_in(2),DSbins(2));
end
% create parameter mesh
[D,sig]=meshgrid(D_in,sig_in);
dDdSig=(D_in(2)-D_in(1))*(sig_in(2)-sig_in(1)); % area element

% compute densities
ln_rho=zeros(size(D,1),size(D,2),W.N);
rho=ln_rho;

   %lnZ=-W.M.ng.*log(W.M.cg)-W.M.na.*log(W.M.ca)...
   %    +gammaln(W.M.ng)+gammaln(W.M.na); % normalization constant
   dt=W.param.timestep;
   b=W.param.blur_beta;
for n=1:W.N
    na=W.M.na(n);
    ng=W.M.ng(n);
    pa=W.M.ca(n)./(2*b*dt*D+sig.^2);
    pg=W.M.cg(n)./(2*dt*D);
    
    
    ln_rho(:,:,n)=...
            log(4*dt*sig/W.M.ca(n)/W.M.ca(n))...
            -gammaln(na)-gammaln(ng)...
            +(na+1)*log(pa)-pa...
            +(ng+1)*log(pg)-pg;        
    %ln_rho(:,:,n)=ln_rho(:,:,n)-max(max(ln_rho(:,:,n)));
    rho(:,:,n)=exp(ln_rho(:,:,n));
end



