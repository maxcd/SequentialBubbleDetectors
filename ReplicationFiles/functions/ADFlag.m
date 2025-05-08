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