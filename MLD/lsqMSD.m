% lsqMSD -- Mean-square-displacement (MSD)-based estimate of diffusion
% constant with motional blur and localization errors.
%
% [Dmsd,S2msd]=lsqMSD(dat,dt,R,NP)
%
% Fit diffusion constant and localization errors to the first NP points of
% an MSD curve. The optimal number of points is usually low, and can be
% estimated by NP ~ floor (  2+2.3*x^0.52 ), where x = sigma^2/D/dt is
% a normalized localization error [1]. Confidence intervals are not
% computed.
%
% WARNING: This file is only intended for use as a cautinary example,
% becasue WHEN PRECISION IS AN ISSUE, THE MSD IS A POOR TO MISERABLE
% ESTIMATOR [2]. A better and equally simple covariance-based estimate has
% been proposed and analyzed by Vestergaard et al [2].
% 
% Input:
% dat   Input data structure, as computed by e.g.,
%       ML1_preprocess_mixed_columns.
% dt    sampling timestep
% R     Blur factor. For images aquired with uniform intensity
%       during a time tE<=dt, R = tE/dt/6 [1].
%
% Output:
% Dmsd  Estimated diffusion constant.
% S2msd Estimated localization error.
% 
%
%
% 1. Xavier Michalet Phys. Rev. E 82, 041914 (2010)
% 2. Vestergaard, Blainey and Flyvbjerg, PRE 89, 022726 (2014)
% 
% Martin Lindén, bmelinden@gmail.com, 2015-04-01

% start of actual code
function [Dmsd,S2msd]=lsqMSD(dat,dt,R,NP)


% compute step length statistics
MSD =zeros(1,NP);
nMSD=zeros(1,NP);
dim=size(dat.x,2); % we will treat each dimension as an independent trajectory
for np=1:NP % loop over step-length interval
    T=dat.end-dat.one+1; % do not assume that dat.T is consistent
    indT=find(T>=np+1);  % at least np+1 positions needed for np steps
    
    for k=1:length(indT)
        rows=dat.one(indT(k)):dat.end(indT(k));
        % sum of square displacements
        MSD(np) = MSD(np)+sum(sum((dat.x(rows((1+np):end),:)-dat.x(rows(1:end-np),:)).^2));
        % sum up number of terms        
        nMSD(np)=nMSD(np)+(length(rows)-np)*dim;
    end
end
MSD=MSD./nMSD;

% < MSD(n) > = 2*D*dt*n + 2*(sigma^2-2*D*dt*R)
P=polyfit(1:NP,MSD,1); % y = P(1)*n + P(2)

Dmsd=P(1)/2/dt;
S2msd= P(2)/2+2*Dmsd*dt*R;






