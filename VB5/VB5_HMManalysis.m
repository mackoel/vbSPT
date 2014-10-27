function res=VB5_HMManalysis(runinputfile)
% res=VB5_HMManalysis(runinputfile)
%
% Run the variational HMM analysis specified in the runinputfile, which
% should be in the current directory. 
% It is also possible to use an options structure, e.g., from
% opt=VB5_getOptions(runinputfile) instead, which might be useful for
% scripting things like parameter sweeps from a single runinput file.
%
% res : structure containing the result of the analysis. The same
% information is also written to the outputfile specified in the runinput
% file or options structure. 

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB5_HMManalysis, runs data analysis in the vbSPT package
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
tstart=tic;
%% read analysis parameters
% if an existing file, generate options structure
if(ischar(runinputfile))
    if(exist(runinputfile,'file')==2)
        opt=VB5_getOptions(runinputfile);
        %disp(['Read runinput file ' runinputfile])
    else
       error(['Cannot find runinput file ' runinputfile '. Please check name and path.'])
    end
elseif(isstruct(runinputfile))
    opt=runinputfile;
    runinputfile=opt.runinputfile; %#ok
    %disp(['Read options structure based on runinput file ' runinputfile ])
else
    error('Could not find runinputfile or interpret runinputfile argument.')
end

% add .mat extension to output file if not present
[outpath,outfile]=fileparts(opt.outputfile);
opt.outputfile=fullfile(outpath,[outfile '.mat']);
% make sure the out folder exist
if(~exist(outpath,'dir'))
   disp(['Creating output folder ' outpath])
   mkdir(outpath)
end
clear outfile outpath;

% construct log file name from outputfile, with extension .log
[logpath,logfile]=fileparts(opt.outputfile);
opt.logfile=fullfile(logpath,[logfile '.log']);
clear logfile logpath;

% start the diary, and clear old entries (the output is overwritten too). 
if(exist(opt.logfile,'file'))
    delete(opt.logfile)
end
diary(opt.logfile);
diary on
VB5_license('VB5_HMManalysis')
disp('----------')
disp([ datestr(now) ' : Starting greedy optimization to find best model.'])
disp(['jobID        : ' opt.jobID])
disp(['runinput file: ' opt.runinputfile])
disp(['input  file  : ' opt.inputfile])
disp(['trj field    : ' opt.trajectoryfield])
disp(['output file  : ' opt.outputfile])
disp(['log file     : ' opt.logfile])
disp('----------')
%% start of actual analysis code
%% set up model structure: prior distributions
%% converge models
% load data
dat=VB5_preprocess(opt);
maxHidden=opt.maxHidden;
timestep=opt.timestep;

%nTot=sum(dat.T-1); % total number of time steps in data
%Ntrj=length(dat.T); % number of trajectories
Witer  =cell(1,opt.runs); % save all models generated in each run

% setup parallel computing
if(opt.parallelize_config)
    if(~isempty(gcp('nocreate'))) % disable existing parpool
        delete(gcp('nocreate'));
    end
    eval(opt.parallel_start)
end

parfor iter=1:opt.runs    
%for iter=1:opt.runs    %%% debug without parfor
    od=0; tx0=0; %#ok
    % Greedy search strategy is probably more efficient than to start over
    % at each model size. We simply start with a large model, and
    % systematically remove the least occupied statate until things start
    % to get worse.
    titer=tic;
    w=VB5_createPrior(opt,maxHidden);
    % initial parameter guess
    Ddt=opt.init_D*timestep;
    Ddt=exp(log(Ddt(1))+diff(log(Ddt))*rand(1,maxHidden)); % log(Ddt) has flat distribution
    
    locErr=opt.init_locErr;
    locErr=exp(log(locErr(1))+diff(log(locErr))*rand(1,maxHidden)); % log(locErr) has flat distribution
    
    ntD=opt.init_tD/timestep;                               % init dwell time interval
    ntD=exp(log(ntD(1))+diff(log(ntD))*rand(1,maxHidden));  % log(ntD) has flat distribution
    A0=diag(1-1./ntD)*eye(maxHidden)+diag(1./ntD/(maxHidden-1))*(ones(maxHidden,maxHidden)-eye(maxHidden));
    
    p0=ones(1,maxHidden)/maxHidden;
    strength=1e4;
    
    w=VB5_newParameters(Ddt,locErr,A0,p0,strength,w);
    %clear Ddt locErr ntD A0 p0 strength
    
    % initial trajectory guess
    w=VB5_trjInit_dat(dat,w,0,0);
    
    % converge largest model
    w0=w;
    w=opt.VBEMfunction(w0,dat,'outputLevel',0,'maxIter',opt.maxIter,...
        'relTolF',opt.relTolF,'tolPar',opt.tolPar);
    %% greedy search
    Witer{iter}{1}=w;
    oneMoreTry=true;
    move=0;
    while(oneMoreTry) % try successive removal of low-occupancy states
        foundImprovement=false;
        move=move+1;
        w0=w; % reference state
        [~,h]=sort(w0.est.Ptot); % order in increasing occupancy
        for k=1:1%length(h) %%% k=1 only is where one decides to only try the least occupied state
            if(w0.N>1)
                %% try to remove a looping state
                w=VB5_removeState(w0,h(k),opt);
                try
                    w=opt.VBEMfunction(w,dat,'outputLevel',0,'maxIter',...
                        opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar);
                catch me
                    disp('VB5_HMManalysis encountered an error:')
                    disp(me.message)
                    disp(me.stack)
                    %%%errFile=['errorlog_VB5_HMManalysis' num2str(rand) '.mat'];
                    %%%save(errFile); % really would like to save this file 
                    rethrow(me);
                end
                Witer{iter}{end+1}=w; % save this attempt
                if(w.F>w0.F) % then this helped, and we should go on
                    w0=w;
                    foundImprovement=true; %#ok
                    %disp(['Iter ' int2str(iter) ': removing state ' int2str(h(k)) ' helped, new size N = ' int2str(w0.N)])
                    break % do not make further attempts at this model size
                end
                % if this did not help, try adding more transition counts and reconverge.
                
                if(w.N>1)
                    tx0=tic;
                    %disp(['Iter ' int2str(iter) ': simple removal did not help. Trying to add some extra transitions'])
                    w=VB5_removeState(w0,h(k),opt);
                    od=max(max(w.M.wB-w.PM.wB)); % largest off-diagonal element
                    w.M.wB=od*(1-eye(w.N));
                    if(isfield(w,'E')) % make sure that the new M field is used
                        w=rmfield(w,'E');
                    end
                    try
                        w=opt.VBEMfunction(w,dat,'outputLevel',0,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar);
                    catch me
                        disp('VB5_HMManalysis encountered an error:')
                        disp(me.message)
                        disp(me.stack)
                        %%%errFile=['errorlog_VB5_HMManalysis' num2str(rand) '.mat'];
                        %%%save(errFile);
                        rethrow(me);
                    end
                    Witer{iter}{end+1}=w; % save this attempt too
                    
                    if(w.F>w0.F)
                        disp(['Iter ' int2str(iter) '. Removing state ' int2str(h(k)) ' (of ' int2str(w0.N) ...
                            ' ) + adding  ' int2str(od) ' extra transitions helped, dF/|F| = ' ...
                            num2str((w.F-w0.F)/abs(w.F)) ', time : ' num2str(toc(tx0)) ' s.']);
                        w0=w;
                        foundImprovement=true; %#ok
                        break % do not make further attempts at this model size
                    else
                        disp(['Iter ' int2str(iter) '. Removing state ' int2str(h(k)) ' (of ' int2str(w0.N) ...
                            ' ) + adding ' int2str(od) ' extra transitions did not help, dF/|F| = ' ...
                            num2str((w.F-w0.F)/abs(w.F)) ', time : ' num2str(toc(tx0)) ' s.']);
                    end
                end
            end
            % give up if W0 could not be improved upon
            if(~foundImprovement)
                oneMoreTry=false;
            end
            %% last attempt: add transition pseudocounts to w0 and reconverge.
            if(w0.N>1)
                tx0=tic;
                w=w0;
                od=max(max(w.M.wB-w.PM.wB)); % largest off-diagonal element
                w.M.wB=od*(1-eye(w.N));
                if(isfield(w,'E')) % make sure that the new M field is used
                    w=rmfield(w,'E');
                end
                try
                    w=opt.VBEMfunction(w,dat,'outputLevel',0,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar);
                catch me
                    disp('VB5_HMManalysis encountered an error:')
                    disp(me.message)
                    disp(me.stack)
                    %%%errFile=['errorlog_VB5_HMManalysis' num2str(rand) '.mat'];
                    %%%save(errFile);
                    rethrow(me);
                end
                Witer{iter}{end+1}=w; % save this attempt too
                if(w.F>w0.F)
                    disp(['Iter ' int2str(iter) '. Adding ' int2str(od) ...
                        ' extra transitions to final model helped : dF/|F| = ' ...
                        num2str((w.F-w0.F)/abs(w.F)) ', time : ' num2str(toc(tx0)) ' s.']);
                else
                    disp(['Iter ' int2str(iter) '. Adding ' int2str(od) ...
                        ' extra transitions to final model did not help : dF/|F| = ' ...
                        num2str((w.F-w0.F)/abs(w.F)) ', time : ' num2str(toc(tx0)) ' s.']);
                end
            end
        end
    end
    
    disp(['Iter ' int2str(iter) '. Finished greedy search in '  num2str(toc(titer)) ' s, with ' int2str(w0.N) ' states.'] )
end
%% collect best models for all sizes
INF=[];
Wbest.F=-inf;
WbestN=cell(1,maxHidden);
bestIter=0;
for k=1:maxHidden
    WbestN{k}.F=-inf;
end
for iter=1:opt.runs
    for k=1:length(Witer{iter})
        w=VB5_sortModel(Witer{iter}{k});
        INF(end+1,1:3)=[iter w.N w.F]; %#ok
        if(w.F>Wbest.F)
            Wbest=w;
        end
        if(w.F>WbestN{w.N}.F)
            WbestN{w.N}=w;
            bestIter=iter;
        end
    end
end
% regenerate unsorted fields
disp('Sorting and regenerating best model and estimates')
if(opt.stateEstimate)
    Wbest=opt.VBEMfunction(Wbest,dat,'estimate','maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar,'outputLevel',0);
else
    Wbest=opt.VBEMfunction(Wbest,dat,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar,'outputLevel',0);
end
parfor k=1:length(WbestN)
    if(WbestN{k}.F>-inf);
        WbestN{k}=opt.VBEMfunction(WbestN{k},dat,'maxIter',opt.maxIter,'relTolF',opt.relTolF,'tolPar',opt.tolPar); %#ok to broadcast opt
    end
end
%
disp(['Best model size: ' int2str(Wbest.N) ', from iteration ' int2str(bestIter) '.'])
%% write results to outputfile
res=struct;
res.options=opt;
res.Wbest=Wbest;
res.WbestN=WbestN;
res.INF=INF;
for k=1:length(WbestN)
    res.dF(k)=WbestN{k}.F-Wbest.F;
end

% saving the models prior to bootstrapping them
disp(['Saving ' opt.outputfile ' after ' num2str(toc(tstart)/60) ' min.']);
save(opt.outputfile,'-struct','res');


%% bootstrapping

if(opt.bootstrapNum>0)
bootstrap = VB5_bsResult(opt, 'HMM_analysis');
res.bootstrap=bootstrap;

% save again after bootstrapping
save(opt.outputfile,'-struct','res');
end

% End parallel computing
if(opt.parallelize_config)
    eval(opt.parallel_end)
end

disp([datestr(now) ' : Finished ' opt.runinputfile '. Total run time ' num2str(toc(tstart)/60) ' min.'])
diary off
end

