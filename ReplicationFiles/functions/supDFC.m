function [reject, pval, tstat, index]=supDFC(y,varargin)
%--------------------------------------------------------------
% PURPOSE: performs supDFC test as in Homm and Breitung (2012)
%--------------------------------------------------------------
% USAGE: [tstat,index] = supDFC(y,d)
% where: y  = a time series vector of dimensions (T,1)
%        d  = optional scalar indicating whether y has been detrended
%             (d=1 in the case with detrending)
%--------------------------------------------------------------
% RETURNS:
%        tstat = scalar value of the test-statistic
%        index = scalar estimate of the break-point
%--------------------------------------------------------------

if nargin > 2
    error('the function supDFC takes at most 2 input arguments')
end

T=length(y);
if nargin ==1 ||  varargin{1}==0
    crit=[1.5762; 1.9327; 2.2685];
elseif varargin{1}==1
    x=[ones(T,1) (1:T)'];
    y=y-x*(x\y);
    crit=[0.9436; 1.3379; 2.0741];
else
    error('second function argument must be logical')
end

r0=floor(.1*T);
dy=diff(y);
test=zeros(T-r0,1);
x=y(1:T-1);
b=x\dy;
e=dy-x*b;
var=(e'*e/(T-2))/(x'*x);
test(1)=b/sqrt(var);
for t0=2:T-r0;
    x0=zeros(t0-1,1);
    x1=y(t0:T-1);
    x=[x0;x1];
    b=x\dy;
    e=dy-x*b;
    var=(e'*e/(T-2))/(x'*x);
    test(t0)=b/sqrt(var);
end;
[tstat, index] =max(test);
    
%if tstat > crit(3)
    %display('H0 is rejected at 1%-level')
%    reject=1;
%    pval=0.01;
if tstat > crit(2)
    %display('H0 is rejected at 5%-level')
     reject=1;
     pval=0.05;
%elseif tstat > crit(1)
%     reject=1;
%     pval=0.1;
else
    %display('H0 is not rejected at 10%-level')
    reject=0; pval=">0.1";
end

%fprintf('The test statistic equals \n')
end