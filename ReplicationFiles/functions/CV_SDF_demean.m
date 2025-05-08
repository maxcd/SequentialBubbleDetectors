%**************************************************************************
%   "Testing for Multiple Bubbles" by Phillips, Shi and Yu (2015)
    
%   In this program, we calculate critical values for the sup 
%   ADF statistic.
% *************************************************************************
 
function [cv_sadf,cv_badf]=CV_SDF_demean(T,swindow0,qe)

m=10000;
dim=T-swindow0+1;

%% %%%% DATA GENERATING PROCESS %%%%%%
SI=1;
randn('seed',SI);   
e=randn(T,m); 
a=1;%T^(-1);  % constant
b =1;% 1/10;  % trend
y=cumsum(e)+a+(1:T)'*b;

%% THE SUP ADF TEST %%%%%%
x=[ones(T,1) (1:T)'];
%x=ones(T,1);% (1:T)'];

badfs=zeros(m,dim); 
%sadf=ones(m,1);
for j=1:1:m
    badf = nan(T-swindow0+1,1);

    % demean
    yy = y(:,j);
    % yy = yy-x*(x\yy);
    for i=swindow0:1:T
      
      % demean in every test window
      yi = yy(1:i)-x(1:i,:)*(x(1:i,:)\yy(1:i));
      badf(i-swindow0+1,1) = ADFlag(yi,0,0);

    end

    badfs(j,:) = badf';
end
sadf=max(badfs,[],2);

cv_sadf=quantile(sadf,qe)';
cv_badf=quantile(badfs,qe)';

%%
% dates = nan(m,1);
% for j=1:1:m
%     if sum(badfs(j,:)' > cv_badf) >= 1
%         ind = find(badfs(j,:)' > cv_badf);
%         dates(j) = ind(1)+swindow0;
%     end
% end
end


