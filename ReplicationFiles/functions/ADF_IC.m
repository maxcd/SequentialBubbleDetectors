function ADFlag = ADF_IC(y, adflag, det, IC)
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
    beta0 = (x1' * x1) \ (x1' * dy);
    eps0 = dy - x1 * beta0;
    T0 = length(eps0);
    npdf0 = sum(-1/2 * log(2 * pi) - 1/2 * (eps0.^2));
    if (IC == 1)  % AIC
        IC0 = -2 * npdf0 / T0 + 2 * length(beta0) / T0;
    elseif (IC == 2) % BIC
        IC0 = -2 * npdf0 / T0 + length(beta0) * log(T0) / T0;
    end
    if (det == 1)
        se0 = (eps0' * eps0) / (t1 - 2);
    end
    if (det == 2)
        se0 = (eps0' * eps0) / (t1 - 3);
    end
    if (det == 0)
        se0 = (eps0' * eps0) / (t1 - 1);
    end
    sig0 = sqrt(diag(se0 * (x1' * x1)^-1));
    ADF0 = beta0(1) / sig0(1);
    IC1 = [];
    ADF1 = [];
    if (adflag > 0)
        IC1 = zeros(adflag, 1);
        ADF1 = zeros(adflag, 1);
        for k = 1:adflag
            t2 = t1 - k;
            xx = x1((k + 1):t1, :);
            dy01 = dy((k + 1):t1);
            x2 = [xx, zeros(t2, k)];
            for j = 1:k
                x2(:, (size(xx, 2) + j)) = dy((k + 1 - j):(t1 - j));
            end
            beta = (x2' * x2) \ (x2' * dy01);
            eps = dy01 - x2 * beta;
            t = length(eps);
            npdf = sum(-1/2 * log(2 * pi) - 1/2 * (eps.^2));
            if (IC == 1)
                IC1(k) = -2 * npdf / t + 2 * size(beta, 1) / t;
            elseif (IC == 2)
                IC1(k) = -2 * npdf / t + size(beta, 1) * log(t) / t;
            end
            if (det == 1)
                se = (eps' * eps) / (t2 - adflag - 2);
            end
            if (det == 2)
                se = (eps' * eps) / (t2 - adflag - 3);
            end
            if (det == 0)
                se = (eps' * eps) / (t2 - adflag - 1);
            end
            sig = sqrt(diag(se * (x2' * x2)^-1));
            ADF1(k) = beta(1) / sig(1);
        end
    end
    ICC = [IC0; IC1];
    ADF = [ADF0; ADF1];
    lag = find(ICC == min(ICC), 1);
    ADFlag = ADF(lag);
    % if (IC == 1)
    %     result = struct('ADF_Statistic_using_AIC', ADFlag);
    % end
    % if (IC == 2)
    %     result = struct('ADF_Statistic_using_BIC', ADFlag);
    % end
end