function [alpha,ln_alpha,lambda_inv,ln_lambda,KL_lv,lambda,nu,sqrtNu]=...
    VB5_ag_vl_expectations(W,showWarning,showAll)
%%
% [alpha,ln_alpha,lambda_inv,ln_lambda,KL_lv,lambda,nu,sqrtNu]=...
%                             VB5_ag_vl_expectations(W,showWarning,showAll)
%
% Numerical evaluation of expectation values of the variational
% distribution q(lambda,v) for the VB5 model.

if(~exist('showWarning','var') || isempty(showWarning))
    showWarning=false;
end
if(~exist('showAll','var') || isempty(showAll))
    showAll=false;
end% global parameters
b=W.param.blur_tau*(1-W.param.blur_tau)-W.param.blur_R;
% tolerance parameters
Nstd_safety=5; % number of estimated standard deviations to include in integrals
relTol=1e-10;    % relative integration tolerance
%%
for state=1:W.N
    showDistribution=showAll;
    %% translate parameters
    n0l=W.PM.nl(state);
    c0l=W.PM.cl(state);
    nl=W.M.nl(state);
    cl=W.M.cl(state);
    n0v=W.PM.nv(state);
    c0v=W.PM.cv(state);
    na=W.M.na(state);
    ca=W.M.ca(state);

    % log density function and derivatives
    logQ=@(l,v)(-nl.*log(l)-cl./l-na.*log(v+b*l)-ca./(v+b*l)...
        -(n0v+1).*log(v)-c0v./v);
    dlogQdl=@(l,v)(  -nl./l     +cl./l.^2  -b*na./(v+b*l)       +b*ca./(v+b*l).^2);
    d2logQdl2 =@(l,v)(nl./l.^2-2*cl./l.^3+b^2*na./(v+b*l).^2-2*b^2*ca./(v+b*l).^3);

    dlogQdv=@(l,v)(  -(n0v+1)./v     +c0v./v.^2-na./(v+b*l)       +ca./(v+b*l).^2);
    d2logQdv2 =@(l,v)((n0v+1)./v.^2-2*c0v./v.^3+na./(v+b*l).^2  -2*ca./(v+b*l).^3);
    d2logQdldv=@(l,v)(                        b*na./(v+b*l).^2-2*b*ca./(v+b*l).^3);    
    %% find maxima a integration domain
    % global maximum in the positive quadrant
    lv0=[cl/(nl+1) c0v/(n0v+1)];
    lv1=fminsearch(@(xx)(-logQ(abs(xx(1)),abs(xx(2)))),lv0);
    lv1=abs(lv1);
    l_max=lv1(1);
    v_max=lv1(2);
    logQref=logQ(l_max,v_max);
    
    % Gaussian approximation of integration boundaries
    hessianMatrix=[d2logQdl2(l_max,v_max) d2logQdldv(l_max,v_max);
                  d2logQdldv(l_max,v_max)  d2logQdv2(l_max,v_max)];
    covarianceMatrix=inv(-hessianMatrix);
    
    stdL=sqrt(covarianceMatrix(1,1));
    stdV=sqrt(covarianceMatrix(2,2));

    l0G=max(0,l_max-12*stdL);
    l1G=l_max+10*stdL;
    v0G=max(0,v_max-12*stdV);
    v1G=v_max+10*stdV;
    %% check for signs of bad scaling        
    corrLV=covarianceMatrix(1,2)/stdL/stdV;    
    if(sqrt(cond(covarianceMatrix))>1000)
        warning(['possible bad scaling: s_max/s_min = '...
            num2str(sqrt(cond(covarianceMatrix))) ...
            ', corrcoeff = ' num2str(corrLV) ])
        showDistribution=showWarning;
    end
    %% unscaled integration     
    tic
    % normalization constant
    zFun=@(ll,vv)(exp(logQ(ll,vv)-logQref));
    Zlv=integral2(zFun,l0G,l1G,v0G,v1G,'relTol',1e-9);
    lnZ=log(Zlv)+logQref; % full normalization constant for exp(logQ).
    
    % <alpha>
    funAlpha=@(ll,vv)(1./(vv+b*ll));
    aFun=@(ll,vv)(funAlpha(ll,vv).*zFun(ll,vv));
    alpha(state)=integral2(aFun,l0G,l1G,v0G,v1G,'relTol',1e-9)/Zlv;

    % <ln alpha>
    lnaFun=@(ll,vv)(log(funAlpha(ll,vv)).*zFun(ll,vv));
    ln_alpha(state)=integral2(lnaFun,l0G,l1G,v0G,v1G,'relTol',1e-9)/Zlv;

    % <gamma> = <1/lambda>
    gFun=@(ll,vv)(1./ll.*zFun(ll,vv));
    lambda_inv(state)=integral2(gFun,l0G,l1G,v0G,v1G,'relTol',1e-9)/Zlv;

    % <ln lambda>
    lngFun=@(ll,vv)(log(ll).*zFun(ll,vv));
    ln_lambda(state)=integral2(lngFun,l0G,l1G,v0G,v1G,'relTol',1e-9)/Zlv;

    if(nargout>=5)
        % <2*D*dt>=<lambda>
        lambdaFun=@(ll,vv)(ll.*zFun(ll,vv));
        lambda(state)=integral2(lambdaFun,l0G,l1G,v0G,v1G,'relTol',1e-9)/Zlv;
    end
    if(nargout>=6)
        % <v> = <sigma^2>
        vFun=@(ll,vv)(vv.*zFun(ll,vv));
        nu(state)=integral2(vFun,l0G,l1G,v0G,v1G,'relTol',1e-9)/Zlv;
    end
    if(nargout>=7)
        % <sig>=<sqrt(v)>
        sigFun=@(ll,vv)(sqrt(vv).*zFun(ll,vv));
        sqrtNu(state)=integral2(sigFun,l0G,l1G,v0G,v1G,'relTol',1e-9)/Zlv;
    end
    tInt=toc;
    % Kullback-Leibler terms
    KL_lv(state)=-lnZ+gammaln(n0l)+gammaln(n0v)-n0l*log(c0l)-n0v*log(c0v)...
            -(nl-n0l)*ln_lambda(state)-(cl-c0l)*lambda_inv(state)...
            +na*ln_alpha-ca*alpha;
    %% plot log density in lambda-v plane
    if(showDistribution)
        ll=linspace(0.9*l0G,1.1*l1G,1e3);
        vv=linspace(0.0*v0G,1.1*v1G,1e3);
        logQref=logQ(l_max,v_max);
        
        [LL,VV]=meshgrid(ll,vv);
        
        figure(1)
        clf
        hold on
        
        %contourf(LL,VV,exp(logQ(LL,VV)-logQref),logspace(-3,0,16),'edgecolor','none')
        contourf(LL,VV,exp(logQ(LL,VV)-logQref),'edgecolor','none')
        plot(l_max,v_max,'k*')
        
        lBox=[l0G l1G l1G l0G l0G];
        vBox=[v0G v0G v1G v1G v0G];
        plot(lBox,vBox,'r')
        
        title(['q(\lambda,\nu), state ' int2str(state) ])
        xlabel('\lambda')
        ylabel('\nu')
        colorbar
        box on
        grid on
        axis([0.95*l0G 1.05*l1G 0.95*v0G 1.05*v1G])
        pause
    end
end

return
%% below: old code for rescaling the integration domain to make
% numericalintegration easier. Leave unused for now.

%% rescale coordinates for integration
%[A,B]=eig(-hessianMatrix);% A'*B*A = - hessianMatrix
%sig=sqrt(inv(B));
%S=sig\A;    % S'*S = -hessianMatrix
S=sqrtm(-hessianMatrix); % S'*S = -hessianMatrix
T=inv(S);

% X = S*(lv-lv_max), lv=lv_max + T*X, dldv = |T|dxdy
funX=@(l,v)(S(1,1)*(l-l_max)+S(1,2)*(v-v_max));
funY=@(l,v)(S(2,1)*(l-l_max)+S(2,2)*(v-v_max));
funL=@(x,y)(l_max+T(1,1)*x+T(1,2)*y);
funV=@(x,y)(v_max+T(2,1)*x+T(2,2)*y);

XX=funX(LL,VV);%S(1,1)*(LL-l_max)+S(1,2)*(VV-v_max);
YY=funY(LL,VV);%S(2,1)*(LL-l_max)+S(2,2)*(VV-v_max);

%xyBoxlv=S*[lBox(1:4)-l_max;vBox(1:4)-v_max];
xyBoxlv=[funX(lBox(1:4),vBox(1:4));funY(lBox(1:4),vBox(1:4))];
xyBox=10*[-1  1  1 -1 ;
    -1 -1  1  1 ];

lvBoxxy=[funL(xyBox(1,:),xyBox(2,:));funV(xyBox(1,:),xyBox(2,:))];


figure(3)
clf
hold on
contourf(XX,YY,exp(logQ(LL,VV)-logQref),logspace(-5,0,16),'edgecolor','none')
plot(xyBoxlv(1,:),xyBoxlv(2,:),'sr-')
plot(xyBox(1,:),xyBox(2,:),'ob-')

% check if xyBox maps to negative lambda/v, in which case more careful
% scaling is needed
if(isempty(find(lvBoxxy(1:end)<0,1))) % then a squre ox is OK!
    xL= xyBox(1,1);
    xU= xyBox(1,2);
    funYL=xyBox(2,1);
    funYU=xyBox(2,3);
    yL=funYL;
    yU=funYU;
    axis equal
    axis([1.1*xL 1.1*xU 1.1*min(yL) 1.1*max(yU)])
    
    figure(1)
    plot(lvBoxxy(1,:),lvBoxxy(2,:),'-bo')
else
    
    % lower bound:
    funY1L=@(x)(xyBox(2,1)*ones(size(x))); % lower limit in xyBox
    funY1U=@(x)(xyBox(2,3)*ones(size(x))); % lower limit in xyBox
    
    [xlv,xilv]=sort(xyBoxlv(1,:));
    ylv=xyBoxlv(2,xilv);
    y3at2=interp1(xlv([1 3]),ylv([1 3]),xlv(2),'linear');
    if(ylv(2)>y3at2)
        funY2L=@(x)(interp1(xlv([1 3]),ylv([1 3]),x,'linear','extrap'));
        funY3L=@(x)(interp1(xlv([3 4]),ylv([3 4]),x,'linear','extrap'));
        funY2U=@(x)(interp1(xlv([1 2]),ylv([1 2]),x,'linear','extrap'));
        funY3U=@(x)(interp1(xlv([2 4]),ylv([2 4]),x,'linear','extrap'));
    else
        funY2L=@(x)(interp1(xlv([1 2]),ylv([1 2]),x,'linear','extrap'));
        funY3L=@(x)(interp1(xlv([2 4]),ylv([2 4]),x,'linear','extrap'));
        funY2U=@(x)(interp1(xlv([1 3]),ylv([1 3]),x,'linear','extrap'));
        funY3U=@(x)(interp1(xlv([3 4]),ylv([3 4]),x,'linear','extrap'));
    end
    funYL=@(x)(max(max(funY1L(x),funY2L(x)),funY3L(x)));
    funYU=@(x)(min(min(funY1U(x),funY2U(x)),funY3U(x)));
    
    %funYL=@(x)(min(max(funY1L(x),funY2L(x)),funY1U(x)));
    %funYU=@(x)(max(min(funY1U(x),funY2U(x)),funY1L(x)));
    
    xL=xyBox(1,1);%min([xyBox(1,:) xyBoxlv(1,:)]);
    if(funYL(xL)>funYU(xL))
        xL=fsolve(@(x)( funYL(x)-funYU(x)),xL);
    end
    xU=xyBox(1,2);%max([xyBox(1,:) xyBoxlv(1,:)]);
    if(funYL(xU)>funYU(xU))
        xU=fsolve(@(x)( funYL(x)-funYU(x)),xU);
    end
    
    xx=linspace(xL,xU,1e4);
    yL=funYL(xx);
    yU=funYU(xx);
    plot(xx,yL,'--k','linew',2)
    plot(xx,yU,'--r','linew',2)
    
    figure(1)
    plot(lvBoxxy(1,:),lvBoxxy(2,:),'ob-')
    lvL=S\[xx;yL]+diag([l_max;v_max])*ones(2,length(xx));
    lvU=S\[xx;yU]+diag([l_max;v_max])*ones(2,length(xx));
    plot(lvL(1,:),lvL(2,:),'--k','linew',2)
    plot(lvU(1,:),lvU(2,:),'--r','linew',2)
    
end

figure(3)
axis equal
axis([1.1*xL 1.1*xU 1.1*min(yL) 1.1*max(yU)])
colorbar
box on
title('ln q(x,y) - rescaled coordinates')
xlabel('x [nm^2]')
ylabel('y [nm^2]')


%% integration 2
dlv=det(T); % area scale factor, cancels in all expectation values

tic
% normalization constant
zFun=@(xx,yy)(exp(logQ(funL(xx,yy),funV(xx,yy))-logQref));
Zxy=integral2(zFun,xL,xU,funYL,funYU,'relTol',1e-9);

% <alpha>
funAlpha=@(xx,yy)(1./(funV(xx,yy)+b*funL(xx,yy)));
aFun=@(xx,yy)(funAlpha(xx,yy).*zFun(xx,yy));
aMean(2)=integral2(aFun,xL,xU,funYL,funYU,'relTol',1e-9)/Zxy;

% <ln alpha>
lnaFun=@(xx,yy)(log(funAlpha(xx,yy)).*zFun(xx,yy));
lnaMean(2)=integral2(lnaFun,xL,xU,funYL,funYU,'relTol',1e-9)/Zxy;

% <gamma>
gFun=@(xx,yy)(1./funL(xx,yy).*zFun(xx,yy));
gMean(2)=integral2(gFun,xL,xU,funYL,funYU,'relTol',1e-9)/Zxy;

% <ln gamma>
lngFun=@(xx,yy)(-log(funL(xx,yy)).*zFun(xx,yy));
lngMean(2)=integral2(lngFun,xL,xU,funYL,funYU,'relTol',1e-9)/Zxy;

% <2*D*dt>=<lambda>
lambdaFun=@(xx,yy)(funL(xx,yy).*zFun(xx,yy));
lambdaMean(2)=integral2(lambdaFun,xL,xU,funYL,funYU,'relTol',1e-9)/Zxy;

% <sig2>=<v>
vFun=@(xx,yy)(funV(xx,yy).*zFun(xx,yy));
vMean(2)=integral2(vFun,xL,xU,funYL,funYU,'relTol',1e-9)/Zxy;

% <sqrt(2*D*dt)>=<sqrt(lambda)>
xStdFun=@(xx,yy)(sqrt(funL(xx,yy)).*zFun(xx,yy));
xStdMean(2)=integral2(xStdFun,xL,xU,funYL,funYU,'relTol',1e-9)/Zxy;

% <sig>=<sqrt(v)>
sigFun=@(xx,yy)(sqrt(funL(xx,yy)).*zFun(xx,yy));
sigMean(2)=integral2(sigFun,xL,xU,funYL,funYU,'relTol',1e-9)/Zxy;

tInt(2)=toc

diff(aMean)/mean(aMean)
diff(lnaMean)/mean(lnaMean)
diff(gMean)/mean(gMean)
diff(lngMean)/mean(lngMean)
diff(lambdaMean)/mean(lambdaMean)
diff(vMean)/mean(vMean)
diff(xStdMean)/mean(xStdMean)
diff(sigMean)/mean(sigMean)

disp('time:')
tInt
tInt/min(tInt)


