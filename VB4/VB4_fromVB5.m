function W=VB4_fromVB5(W)
% function W=VB4_fromVB5(W)
%
% convert a VB5 model to the corresponding VB4 one, using the priors on v
% as priors on alpha in the process.

W.PM.ng=W.PM.nl;
 W.M.ng=W.M.nl;
 try
     W.E.ng=W.E.nl;
 end
 
W.PM.cg=W.PM.cl;
 W.M.cg=W.M.cl;
try
    W.E.cg=W.E.cl;
end

W.PM.na=W.PM.nv;
W.PM.ca=W.PM.cv;

W.PM=rmfield(W.PM,{'nl','cl','nv','cv'});
 W.M=rmfield(W.M,{'nl','cl'});

W.param.blur_beta=W.param.blur_tau*(1-W.param.blur_tau)-W.param.blur_R;
