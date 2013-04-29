%% VB-HMM analysis parameter file generated by vbSPTgui %%

% Fredrik Persson & Martin Linden 2012-06-13

% This is a GUI generated HMM analysis runinput file, which specifies
% everything the code needs to know about how to analyze a particular data
% set.
% To run the HMM analysis manually type:
% >> VB3_HMManalysis('runinputfilename')


% Inputs
inputfile = './InputData/testdata_vbSPT.mat';
trajectoryfield = 'finalTraj';

% Computing strategy
parallelize_config = 1;
parallel_start = 'matlabpool open';  % executed before the parallelizable loop.
parallel_end = 'matlabpool close'; % executed after the parallelizable loop.

% Saving options
outputfile = './Results/testresult_vbSPT_HMM_normal.mat';
jobID = 'a short string to describe this job';

% Data properties
timestep = 0.003;     % in [s]
dim = 1;
trjLmin = 2;

% Convergence and computation alternatives
runs = 25;
maxHidden = 10;

% Evaluate extra estimates including Viterbi paths
stateEstimate = 0;

maxIter = [];    % maximum number of VB iterations ([]: use default values).
relTolF = 1e-8;  % convergence criterion for relative change in likelihood bound.
tolPar = [];     % convergence criterion for M-step parameters (leave non-strict).

% Bootstrapping
bootstrapNum = 100;
fullBootstrap = 1;

% Limits for initial conditions
init_D = [0.1, 10]*1e6;   % interval for diffusion constant initial guess [length^2/time] in same length units as the input data.
init_tD = [2, 20]*timestep;     % interval for mean dwell time initial guess in [s].
% It is recommended to keep the initial tD guesses on the lower end of the expected spectrum.

% Prior distributions
% The priors are generated by taking the geometrical average of the initial guess range.
% units: same time units as timestep, same length unit as input data.
prior_D = 1e6;         % prior diffusion constant [length^2/time] in same length units as the input data.
prior_Dstrength = 5;   % strength of diffusion constant prior, number of pseudocounts (positive).
prior_tD = 10*timestep;      % prior dwell time in [s]. Must be greater than timestep (recommended > 2*timestep)
prior_tDvar = 100*prior_tD;   % variance of prior dwell times [s]. 

