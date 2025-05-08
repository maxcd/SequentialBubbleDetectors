function result = SDF_demean_rec(y, r0, det)

T = length(y);

if nargin < 2
    r0 = 0.01 + 1.8/sqrt(T);
end



    
    
    swindow0 = floor(r0 * T);
    dim = T - swindow0 + 1;
    badfs = zeros(dim, 1);


%    detrend recursively
    dy = y(2:end) - y(1:end-1);
    dy_adj=zeros(T-2,1);        % adjust for a trend recursively
    for j=2:T-1
        dy_adj(j-1)=sqrt((j-1)/j)*(dy(j)-mean(dy(1:j-1)));
    end
    y_adj=[zeros(2, 1); cumsum(dy_adj)];


    for i = swindow0:T
        badfs(i - swindow0 + 1) = ADFlag(y_adj(3:i), 0, 0);
    end
    sadf = max(badfs);

    result = struct('badfs', badfs, 'sadf', sadf);
        %'bsadfs', bsadfs, 'gsadf', gsadf);
end


function ADF1 = ADFlag(y, adflag, det)
    t0 = length(y);
    t1 = length(y) - 1;
    const = ones(t1, 1);
    trend = (1:t1)';
    y1 = y(1:t1);
    dy = y(2:t0) - y(1:t1);
    x = y1;
    if (det == 1)
        x = [x, const];
    end
    if (det == 2)
        x = [x, const, trend];
    end
    if (det == 0)
        x = x;
    end
    x1 = x;

    ADF1 = zeros(length(adflag), 1);
    for k = adflag
        ii = 1;
        t2 = t1 - k;
        xx = x1((k + 1):t1, :);
        dy01 = dy((k + 1):t1);
        x2 = [xx, zeros(t2, k)];
        for j = 1:k
            x2(:, (size(xx, 2) + j)) = dy((k + 1 - j):(t1 - j));
        end
        ncoeff = 1+det+k;
        ixx = (x2' * x2) \ speye(ncoeff );
        beta = ixx * (x2' * dy01);
        eps = dy01 - x2 * beta;
        se = (eps' * eps) / (t2 - ncoeff);
        sig = sqrt(diag(se * ixx));
        ADF1(ii) = beta(1) / sig(1);
        ii = ii +1;
    end
    


end