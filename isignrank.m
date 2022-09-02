function [prob, h, cred_bounds] = isignrank(y,x,varargin)
% ISIGNRANK Bayesian Wilcoxon signed rank sum test based on the Imprecise Dirichlet Process.
%   Prob = ISIGNRANK(Y,X) computes the lower and upper posterior probability
%   of the  hypothesis P(Z>=-Z)>1/2, where Z=Y-X.
%   Prob(1,1) is the lower and Prob(2,1) is the upper.
% 
%   [Prob,H] = ISIGNRANK(...) performs a hypothesis test for P(Z>=-Z)>1/2
%   with posterior probability 0.95.
%   H=1 indicates that the hypothesis P(Z>=-Z)>1/2 is true with posterior
%   probability greater than 0.95, i.e., Prob(1,1)>0.95. H=0 indicates
%   P(Z>=-Z)>1/2 is not true with posterior greater than 0.95, i.e.,
%   Prob(2,1)<0.95. Finally, H=2 indicates an indeterminate instance, i.e.,
%   we cannot decide if  P(Z>=-Z)>1/2 is true with posterior probability
%   greater than 0.95. This means that the posterior inferences are prior
%   dependent, i.e., Prob(1,1)<0.95 and  Prob(2,1)>0.95.
% 
%   [Prob,H] = ISIGNRANK(...,'alpha',ALPHA) returns the result of the hypothesis
%   test with posterior probability gretaer than 1-ALPHA.
% 
%   [Prob,H] = ISIGNRANK(...,'s',sval) sets the value of the prior strength of 
%   the Dirichlet Process s to sval.
%   The default value is s=(sqrt(17)-3)/2. For sval=0 it coincides with the 
%   Bayesian bootstrap.
% 
%   [Prob,H] = ISIGNRANK(...,'method',M) computes the posterior probability if M is
%   'exact', or uses a normal approximation if M is 'approximate'.  If you
%   omit this argument, ISIGNRANK uses the exact method for small samples and
%   the approximate method for larger samples.
% 
%   [Prob,H] = ISIGNRANK(...,'nsamples',N) sets the number of samples generated
%   from the Dirichlet distribution to compute the posterior probabilities.
%   Default value is 200000.
% 
%   [Prob,H] = ISIGNRANK(...,'display','off') does not show the plot.
% 
%   [Prob,H] = ISIGNRANK(...,'rope',val) introduces a (symmetric) Region of Practical
%   Equivalence (ROPE) around 1/2, i.e., [1/2-val,1/2+val].
% 
% 
%   [Prob,H] = ISIGNRANK(...,'tail',TAIL) performs the test specified by
%   TAIL:
%     'right' -- evaluates the hypothesis  P(Z>=-Z)>1/2. This is the
%                Bayesian improved version of ranksum(y,x,'tail','right') (default value).
%     'left'  -- evaluates the hypothesis  P(Y <= X)>1/2. This is the
%                Bayesian improved version of ranksum(y,x,'tail','left').
%     'both'  -- performs a two-sided Bayesian test, i.e., H=1 if 1/2 is
%                not included in the 1-ALPHA lower and upper HPD credible
%                intervals. H=0 if 1/2 is included in the 1-ALPHA lower and
%                upper HPD credible. H=2 otherwise, indeterminate case.
%                This is the Bayesian improved version of signrank(y,x).
% 
%   [Prob] = ISIGNRANK(...,'tail','neighbour','bound',[v1 v2]) computes the
%   integral  between v1 and v2 of the  lower and upper distribution
%   of P(Z>=-Z). Note that, 0<=v1<v2<=1.
% 
%   [Prob, H, Cred_bounds] = ISIGNRANK(...) returns the lower and upper
%   posterior credible interval at level 1-ALPHA, where
%   [cred_bounds(1,1) cred_bounds(1,2)] is the credible interval of the
%   lower distribution of P(Z>=-Z) and [cred_bounds(2,1) cred_bounds(2,2)]
%   is the credible interval of the upper distribution.
% 
%   Examples:
%   x=randn(10,1);
%   y=randn(10,1);
%   [prob,h]=isignrank(y,x)
% 
%   x=randn(10,1);
%   y=randn(10,1)+4;
%   [prob,h]=isignrank(y,x)
% 
%   x=randn(10,1);
%   y=randn(10,1);
%   [prob,h,cred_bounds]=isignrank(y,x,'tail','both')
% 
%   See also ISIGNTEST, IRANKSUM.
% 
%   References:
%      [1] A. Benavoli, F. Mangili, F. Ruggeri and M. Zaffalon
%         "A Bayesian Wilcoxon signed-rank test based on the Dirichlet
%         process" accepted to ICML 2014
%
%
%  This is the matalab implementation of the ISIGNRANK test and
%  is released under the Gnu Public License (GPL). 
%  It is free software and has no warranty. Read the license in COPYING.
% 
%  Copyright (c) 2014 IPG IDSIA
%  alessio@idsia.ch or benavoli@gmail.com
%  IDSIA,  Galleria 2,  6928 Manno,  Switzerland




% addpath('utilities')
prob=[];
%define defaults option
options = struct('alpha',0.05,'s',(sqrt(17)-3)/2,'method','exact','tail','right','nsamples',200000,'display','on','rope',0,'bounds',[0.45 0.55]);

%read the acceptable names
optionNames = fieldnames(options);

if nargin < 2
    error('isignrank:MissingInput', ...
        'Two input vectors are required');
end


if nargin > 2
    % count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
        error('isignrank:BadFormat options needs propertyName/propertyValue pairs')
    end
    
    for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
        inpName = lower(pair{1}); % make case insensitive
        
        if any(strmatch(inpName,optionNames))
            % overwrite options
            options.(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name',inpName)
        end
    end
    
end

method=options.method;
alphav=options.alpha;
tail=options.tail;
s=options.s;
nsamples=options.nsamples;
display=options.display;
rope=options.rope;
VecBound=options.bounds;

z = y(:) - x(:);

% Remove missing data
z(isnan(z)) = [];
if (length(z)==0)
    error('isignrank:noData','No data remaining after removal of NaNs.');
end

nodiff = find(z == 0);
z(nodiff) = [];
 n = length(z);



% Method check
if n>250
    method = 'approximate';
    warning('isignrank:Method',...
        'Data Size is too large. METHOD is set to ''approximate''.')
end

X=repmat(z,1,n);
Y=repmat(-z',n,1);

A = heaviside(X-Y); 
% Au=[[Al ones(nx,1)];[ones(1,ny) 1]];



if strcmp(method,'exact')
    
    data1 = drchrnd([ones(1,n) s],nsamples);
    
    
    data12l=dot(data1(:,1:end-1)*A,data1(:,1:end-1),2);
    
    data12u=data12l+data1(:,end).*(2-data1(:,end));
    
    
    
else
    %Use Normal approximation
    
    CovW = (ones(n,n) + diag(ones(1,n))) /((s + n)*(s + n + 1));
    
    %For the lower
    %compute mean
    meanl=trace(A*CovW);
    
    %compute variance
    den= (s + n + 3)*(s + n + 2)*(s + n + 1)*(s + n);
    mu0=1/den;
    mu1=2/den;
    mu2=4/den;
    mu3=6/den;
    mu4=24/den;
    a=diag(A);
    As=(A+A)/2;
    a4=a'*a; a2=(2*trace(As^2)+trace(As)^2-3*a'*a); a3=4*trace(diag(a)*(As-diag(a))*ones(n,n));
    a1a=2*(trace(diag(a)*ones(n,n)*As*ones(n,n))-2*trace(diag(a)*(As)*ones(n,n))+a'*a)-trace(diag(a)*((diag(a)*ones(n,n))'-diag(a))*ones(n,n))-(trace(As)^2-a'*a);
    a1b=4*trace((As-diag(a))*(As-diag(a))*ones(n,n))-2*(2*trace(As^2)-2*a'*a);
    a1=(a1a+a1b);
    a0=(trace(As*ones(n,n)*As*ones(n,n))-a1-a2-a3-a4);
    Varl=mu4*a4+mu3*a3+mu2*a2+mu1*a1+mu0*a0-meanl^2;
    
    
    
    
    
    %For the Upper
    %compute mean
    meanu=meanl+(s^2 +2*s*n+s)/((s+n)*(s+n+1));
    
    %compute variance
    Q= (ones(n,n) + diag(ones(1,n)))*s*(5+2*n+s)/((n+s)*(1+n+s)*(2+n+s)*(3+n+s));
    Varu=Varl+2*n*(1+n)*s*(3+6*n+2*n^2 +3*s+2*n*s)/((n+s)^2*(1+n+s)^2*(2+n+s)*(3+n+s))+2*trace(A*Q)-2*meanl*(s^2 +2*s*n+s)/((s+n)*(s+n+1));
    
    
    
    
    data12l=normrnd(meanl,sqrt( Varl),nsamples,1);
    data12u=normrnd(meanu,sqrt( Varu),nsamples,1);
    
end

%compute credible interval for the lower
SD = sort(data12l);
xll = SD(ceil(alphav*nsamples/2)); %left bound
xlr = SD(ceil((1-alphav/2)*nsamples)); %right bound

%compute credible interval for the upper
SD = sort(data12u);
xul = SD(ceil(alphav*nsamples/2)); %left bound
xur = SD(ceil((1-alphav/2)*nsamples)); %right bound


if strcmp(tail,'both')
    
    
    if (xlr<0.5-rope && xur<0.5-rope) || (xll>0.5+rope && xul>0.5+rope)
        h=1; %hypotheis H1 is accepted
    elseif  ((xll<0.5-rope && xlr>=0.5+rope) && (xul<0.5-rope && xur>=0.5+rope)) 
        h=0; %hypotheis H0 is accepted
    else
        h=2; %no enough information to discrimnate hypotheses H0  and H1
    end
    
    prob=[];
    
elseif strcmp(tail,'right')
    areal=length(find(data12l>0.5+rope))/nsamples;
    areau=length(find(data12u>0.5+rope))/nsamples;
    
    if areal>1-alphav && areau>1-alphav
        h=1; %hypotheis H1 is accepted
    elseif  (areal<=1-alphav && areau>1-alphav) ||  (areal>1-alphav && areau<=1-alphav)
        h=2; %no enough information to discriminate hypotheses H0  and H1
    else
        h=0;  %hypotheis H0 is accepted
    end
    prob=[areal;areau];
    
elseif strcmp(tail,'left')
    areau=length(find((1-data12l)>0.5+rope))/nsamples;
    areal=length(find((1-data12u)>0.5+rope))/nsamples;
    
    
    if areal>1-alphav && areau>1-alphav
        h=1; %hypotheis H1 is accepted
    elseif  (areal<=1-alphav && areau>1-alphav) ||  (areal>1-alphav && areau<=1-alphav)
        h=2; %no enough information to discrimnate hypotheses H0  and H1
    else
        h=0;  %hypotheis H0 is accepted
    end
    prob=[areal;areau];
    
elseif strcmp(tail,'neighbour')
    areal=length(find(data12l>VecBound(1) & data12l<VecBound(2)))/nsamples;
    areau=length(find(data12u>VecBound(1) & data12u<VecBound(2)))/nsamples;
    
    
    if areal>1-alphav && areau>1-alphav
        h=0; %hypotheis H0 is accepted
    elseif  (areal<=1-alphav && areau>1-alphav) ||  (areal>1-alphav && areau<=1-alphav)
        h=2; %no enough information to discrimnate hypotheses H0  and H1
    else
        h=1;  %hypotheis H1 is accepted
    end
    prob=[areal;areau];
    
end

cred_bounds=[xll, xlr;xul xur];


if strcmp(display,'on')
    
    nylabel='Probability';
    nxlabel=('P(Z > -Z'')');
    
    
    if strcmp(tail,'both')
        
        fh=figure;
        hpl=plot_twoside(data12l,xll,xlr,fh,'k',1,nylabel,nxlabel);
        hold on
        hpu=plot_twoside(data12u,xul,xur,fh,'r',2,nylabel,nxlabel);
        legend([hpl hpu],'lower','upper')
        
        
    elseif strcmp(tail,'right')
        
        fh=figure;
        xlim([0,1])
        hold on
        hpl=plot_oneside(data12l,areal,1,fh,'AreaLower=','k',0.15,nylabel,nxlabel);
        hold on
        hpu=plot_oneside(data12u,areau,1,fh,'AreaUpper=','r',0.75,nylabel,nxlabel);
        legend([hpl hpu],'lower','upper')
        
    elseif strcmp(tail,'left')
        nxlabel=('P(Z < -Z'')');
        fh=figure;
        xlim([0,1])
        hold on
        hpl=plot_oneside(data12l,areau,-1,fh,'AreaUpper=','k',0.15,nylabel,nxlabel);
        hold on
        hpu=plot_oneside(data12u,areal,-1,fh,'AreaLower=','r',0.75,nylabel,nxlabel);
        legend([hpl hpu],'upper','lower')
        
    elseif strcmp(tail,'neighbour')
        
        fh=figure;
        xlim([0,1])
        hold on
        hpl=plot_around(data12l,VecBound(1),VecBound(2),fh,'k',0.15,nylabel,nxlabel,areal,'AreaLower=');
        hold on
        hpu=plot_around(data12u,VecBound(1),VecBound(2),fh,'r',0.75,nylabel,nxlabel,areau,'AreaUpper=');
        legend([hpl hpu],'lower','upper')
    end
    
end


end



