%**************************************************************************
%   "Testing for Multiple Bubbles" by Phillips, Shi and Yu (2015)
    
%   In this program, we calculate critical values for the sup 
%   ADF statistic.
% *************************************************************************
 
function [cv_sadf,cv_badf]=CV_SADF(T,swindow0,qe)

m=10000;
dim=T-swindow0+1;

%% %%%% DATA GENERATING PROCESS %%%%%%
SI=1;
randn('seed',SI);   
e=randn(T,m); 
a=T^(-1);
y=cumsum(e+a);

%% THE SUP ADF TEST %%%%%%

badfs=zeros(m,dim); 
%sadf=ones(m,1);
parfor j=1:1:m
    badf = nan(T-swindow0+1,1);
    for i=swindow0:1:T
      %badfs(j,i-swindow0+1)= ADFlag(y(1:i,j),0,1);%ADF_FL(y(1:i,j),0,1);
      badf(i-swindow0+1,1) = ADFlag(y(1:i,j),0,1);
    end
    badfs(j,:) = badf';
end
sadf=max(badfs,[],2);

cv_sadf=quantile(sadf,qe)';
cv_badf=quantile(badfs,qe)';

% %%
% dates = nan(m,1);
% for j=1:1:m
%     if sum(badfs(j,:)' > cv_badf) >= 1
%         ind = find(badfs(j,:)' > cv_badf);
%         dates(j) = ind(1)+swindow0;
%     end
% end
end


