function [alpha,ln_alpha,lambda_inv,ln_lambda,KL_lv,lambda,v,sqrtv]=...
    VB5_ag_vl_expectations_beta0(W)
%%
% [alpha,ln_alpha,lambda_inv,ln_lambda,KL_lv,lambda,nu,sqrtv]=...
%                             VB5_ag_vl_expectations_bet0(W)
%
% Numerical evaluation of expectation values of the variational
% distribution q(lambda,v) for the VB5 model.

% global parameters
b=W.param.blur_tau*(1-W.param.blur_tau)-W.param.blur_R;
if(b>0)
    warning([' beta = ' num2str(beta) ...
        '. VB5_ag_vl_expectations_beta0 only applies to the case beta=0.'])
end
%%
for state=1:W.N    
    %% translate parameters
    n0l=W.PM.nl(state);
    c0l=W.PM.cl(state);
    nl=W.M.nl(state);
    cl=W.M.cl(state);
    n0v=W.PM.nv(state);
    c0v=W.PM.cv(state);
    na=W.M.na(state);
    ca=W.M.ca(state);

    nv=n0v+na;
    cv=c0v+ca;
    
    %% exact integration     
    % normalization constant
    lnZ(state)=gammaln(nl)-nl*log(cl)+gammaln(nv)-(nv)*log(cv);
    
    % <alpha>
    alpha(state)=nv/cv;
    v_inv(state)=alpha(state);
    
    % <ln alpha> = - <ln v>
    ln_alpha(state)=psi(nv)-log(cv);
    ln_v(state)=-ln_alpha(state);
    
    % <gamma> = <1/lambda>
    lambda_inv(state)=nl/cl;

    % <ln lambda>
    ln_lambda(state)=-(psi(nl)-log(cl));

    % <2*D*dt>=<lambda>
    lambda(state)=cl/(nl-1);
    
    % <v> = <sigma^2>
    v(state)=cv/(nv-1);
    
    % <sig>=<sqrt(v)>
    sqrtv(state)=0.5*log(pi*cv)+gammaln(nv-0.5)-(nv-1)*log(4)-2*gammaln(nv);
    
    % Kullback-Leibler terms
    KL_lv(state)=-lnZ+gammaln(n0l)+gammaln(n0v)-n0l*log(c0l)-n0v*log(c0v)...
            -(nl-n0l)*ln_lambda(state)-(cl-c0l)*lambda_inv(state)...
            -(nv-n0v)*ln_v(state)     -(cv-c0v)*v_inv(state);
%    KL_lv3(state)=-gammaln(nl)+gammaln(n0l)+n0l*log(cl/c0l)...
%                  -gammaln(nv)+gammaln(n0v)+n0v*log(cv/c0v)...
%                  +(nl-n0l)*(psi(nl))-(cl-c0l)*nl/cl...
%                  +(nv-n0v)*(psi(nv))-(cv-c0v)*nv/cv;
%    KL_lv2(state)=+n0l*(log(cl/c0l))-gammaln(nl)+gammaln(n0l)+(nl-n0l)*psi(nl)-nl*(1-c0l/cl)...
%                  +n0v*(log(cv/c0v))-gammaln(nv)+gammaln(n0v)+(nv-n0v)*psi(nv)-nv*(1-c0v/cv);
end
%disp('beta0 KL differences:')
%disp(num2str(KL_lv-KL_lv2))
%disp(num2str(KL_lv3-KL_lv2))
%disp(num2str(KL_lv-KL_lv3))
%disp('---------------------')

