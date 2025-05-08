function [reject, tstat, adf, pval]  = supDF(y,varargin)
%--------------------------------------------------------------
% PURPOSE: performs supDF test as in Homm and Breitung (2012)
%--------------------------------------------------------------
% USAGE: tstat = supDF(y,d)
% where: y  = a time series vector of dimensions (T,1)
%        d  = optional scalar indicating whether y has been detrended
%             (d=1 in the case with detrending)
%--------------------------------------------------------------
% RETURNS:
%        tstat = scalar value of the test-statistic
%---------------------------------------------------------------

if nargin > 2
    error('the function supDF takes at most 2 input arguments')
end

T=length(y);
if nargin ==1 ||  varargin{1}==0
    crit=[2.4152; 2.7273; 3.3457];
elseif varargin{1}==1
    x=[ones(T,1) (1:T)'];
    y=y-x*(x\y);
    crit=[0.5921; 0.8726; 1.4176];
else
    error('second function argument must be logical')
end
    
r0=floor(0.1*T);
adf=zeros(T-r0,1);     % Sequential DF-Test a r0
for i=r0+1:T
     z=y(1:i-1);
     dz=y(2:i)-y(1:i-1);
     b=z\dz;
     u=dz-z*b;
     sig2=u'*u/(i-2);  
     adf(i-r0)=b*norm(z)/sqrt(sig2);
end
tstat=max(adf);

%if tstat > crit(3)
%    display('H0 is rejected at 1%-level')
%    reject=1;
%    pval=0.01;
if tstat > crit(2)
    %display('H0 is rejected at 5%-level')
     reject=1;
     pval=0.05;
%elseif tstat > crit(1)
 %    reject=1;
 %    pval=0.1;
    %display('H0 is rejected at 10%-level')
else
    %display('H0 is not rejected at 10%-level')
    reject=0; pval=">0.1";
end

%fprintf('The test statistic equals ')
end