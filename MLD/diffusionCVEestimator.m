% diffusionCVEestimator -- Covariancebased estimation of diffusion constant
% from trajectories with motional blur and localization errors.
%
% [Dvbf,S2vbf]=diffusionCVEestimator(dat,dt,R,sigma2)
%
% Computes the covariance based diffusion constant estimators by
% Vestergaard, Blainey, and Flyvbjerg (VBF) [1] on a multi-trajectory data
% set. Two implementations are provided depending on input.
%
% a) if a localization error sigma2 (sigma^2 in VBF's notation) is given,
% only the diffusion constant is estimated, using eq. (16) in [1].
% b) if sigma square is empty or not given, then both diffusion constant
% and localization error is estimated, using eqs. (14,15).
%
% In both cases, the arithmetic averages over all trajectories are used.
% 
% Input:
% dat   Input data structure, as computed by e.g.,
%       ML1_preprocess_mixed_columns.
% dt    sampling timestep
% R     Blur factor. For images aquired with uniform intensity
%       during a time tE<=dt, R = tE/dt/6 [1].
% sigma2 Optional localization error of the data.
%
% Output
% Dvbf  Estimated diffusion constant.
% S2vbf Estimated localization error (only in case b).
%
% Ref. 1 also describes confidence interval estimates and consistency
% checks for diffusive motion, which are currently not implemented.
%
% 1. Vestergaard, Blainey and Flyvbjerg, PRE 89, 022726 (2014)
% 
% Martin Lindén, bmelinden@gmail.com, 2015-04-01

% start of actual code
function [Dvbf,S2vbf]=diffusionCVEestimator(dat,dt,R,sigma2)

% decide what to estimate
if( exist('sigma2','var') && ~isempty(sigma2))
   sigma2IsKnown=true; 
else
    sigma2IsKnown=false;
end

% compute step length statistics
ind3=find(dat.T>=3); % at least two steps are needed to estimate
nDX2 =dat.dim*sum(dat.T(ind3)-1); % total number of steps
nDXp1=dat.dim*sum(dat.T(ind3)-2); % total number of step pairs

DX2mean=0;      % < dx(t)^2 >
DXtp1mean=0;    % < dx(t)*dx(t+1) >

for k=1:length(ind3)
    rows=dat.one(ind3(k)):dat.end(ind3(k));
    DX=diff(dat.x(rows,dat.xInd),1,1); % step lengths
    
    DX2mean=DX2mean+sum(sum(DX.^2))/nDX2;
    DXtp1mean=DXtp1mean+sum(sum(DX(1:end-1,:).*DX(2:end,:)))/nDXp1;
    
end
% alternative
dx=diff(dat.x(:,dat.xInd),1,1);
dx(dat.end,:)=0;

dx2mean=sum(sum(dx.^2))/nDX2; % dx(dat.end,:) do not contribute to the sum
dxtp1mean=sum(sum(dx(1:end-1,:).*dx(2:end,:)))/nDXp1; 
% check that the numbers are the same!
keyboard

if (sigma2IsKnown)  % then estimate only D 
    Dvbf = (DX2mean-2*sigma2)/(2*(1-2*R)*dt);
    S2vbf=[];
else                % estimate D and sigma2
    Dvbf = DX2mean/2/dt+DXtp1mean/dt;
    S2vbf= R*DX2mean+(2*R-1)*DXtp1mean;
end



