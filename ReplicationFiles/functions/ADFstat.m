function stat = ADFstat(y, p, det)

% p = number of lags in the test equation (in first differences)
p=p;

T = length(y);
N = T-p-1;                      % effective sample size
dy  = y(2:end) - y(1:end-1);    % first difference
X   = zeros(N,p+1);             % p+1 regressors in the test equation
X(:,1) = y(1+p:T-1);            % lagged dependent variable as first regressor
if p>0                          % add p lags of differences
    for pp=p:-1:1
        X(:,pp+1) = dy(1+(p-pp):N+p-pp);
    end
end
% add deterministic components
switch det
    case 0  % nothing
    case 1  % constant
        X = [X, ones(N, 1)];
    case 2  % constant and trend
        X = [X, ones(N, 1), (1:N)'];
end
% dependent variable: first difference
Y = dy(1+p:N+p);

b= (X'*X)\(X'*Y);
e = Y -X*b;
ncoef = length(b);
df = N-ncoef;
SSE = sum(e.^2);
sig2 = SSE/df;
se = sqrt(sig2./sum(X(:,1).^2) );
stat = b(1)/se;

end
