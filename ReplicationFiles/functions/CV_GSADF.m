%**************************************************************************
%   "Testing for Multiple Bubbles" by Phillips, Shi and Yu (2015)
    
%   In this program, we calculate critical values for the generalized sup 
%   ADF statistic.
% *************************************************************************
 
function [cv_gsadf,cv_bsadf]=CV_GSADF(T,swindow0,qe)


m=10000;
dim=T-swindow0+1;

%% %%%% DATA GENERATING PROCESS %%%%%%
SI=1;
randn('seed',SI);   
e=randn(T,m); 
a=T^(-1);
y=cumsum(e+a);

%% THE GENERALIZED SUP ADF TEST %%%%%%

gsadf=zeros(m,1);  
bsadfs=zeros(m,dim); 
parfor j=1:m
    
    max_rwadft = zeros(1, dim);
    for r2=swindow0:1:T
        
        dim0=r2-swindow0+1;
        
        rwadft=zeros(dim0,1);
        for r1=1:1:dim0
            rwadft(r1)= ADFlag(y(r1:r2,j),0,1);%ADF_FL(y(r1:r2,j),0,1);  % two tail 5% significant level
        end
        %sadfs(j,r2-swindow0+1)=max(rwadft);
        max_rwadft(r2-swindow0+1) = max(rwadft);
    end
    bsadfs(j,:) = max_rwadft;
    gsadf(j)    = max(max_rwadft);
end

cv_gsadf=quantile(gsadf,qe)';
cv_bsadf=quantile(bsadfs,qe)';


