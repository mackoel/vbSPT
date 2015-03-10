function dat=ML1_preprocess(cellData,dim)
% dat=ML1_preprocess(cellData,dim)
%
% assemble single particle diffusion data in a form for fast EM iterations
%
% cellData  : a cell vector of trajectories, where each row is 
%             cellData{k}(t,:) = [x1 x2 ...xd e1 e2 ... ed]
%             xi : positions, ei: errors (standard deviations)
% optional parameters with cell vector input (ignored otherwise):
% dim       : specification of data dimension (d) when input is cell
%             vector; it is then assumed that the first dim columns are
%             positions, and the next dim columns are errors. 
%             Default: use all columns.
% output: data dstructure dat, with fields
% dat.x, dat.v     : positions and errors (variances)
% dat.dim   : dimension
% dat.T     : vector or trajectory lengths
% dat.one, dat.end : indices of starts and ends of all trajectories.
% dat.Yone, dat.Yend: indices to the unknown true trajectories start- and end-points
% (hidden trajectories contain one extra position).
%
% ML 2015-03-10
%% parse input

X=cellData;
if(~exist('dim','var'))
    columns=size(X{1},2);
    dim=floor(columns/2);
    warning(['ML1_preprocess: data dimension not specified, using dim = ' int2str(dim)])
end
%% assemble output structure
dat=struct;
dat.dim=dim;

% count output size
T=zeros(size(X));
for k=1:length(X)
    T(k)=size(X{k},1);
end
if(~isempty(find(T<2,1)))
   warning('VB5_preprocess: data contains traces with no steps.')
end
% data stacking: pack a zero-row between every trajectory to match sizes of
% data x(t) and diffusive path y(t).
dat.T=T;
dat.x=zeros(sum(T),dim);
dat.v=zeros(sum(T),dim);
dat.one=zeros(1,length(X),'double');
dat.end =zeros(1,length(X),'double');

ind=1;
for k=1:length(X)
    x =X{k}(:,1:dim);
    dx=X{k}(:,dim+(1:dim));
    Tx=size(x,1);
    dat.one(k)=ind;    
    dat.end(k)  =ind+Tx-1;
    ind=ind+Tx;
    dat.x(dat.one(k):dat.end(k),1:dim)=x;
    dat.v(dat.one(k):dat.end(k),1:dim)=dx.^2;
end

% indices to the hidden trajectory, which has one extra position...
% indices to hidden trajectory start- and end-points
dat.Yone=dat.one+(0:length(dat.one)-1);
dat.Yend=dat.end+(1:length(dat.one)  );



