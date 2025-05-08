function result = GSADF_IC(y, adflag, det, IC, parallel)
    t = length(y);
    r0 = 0.01 + 1.8/sqrt(t);
    swindow0 = floor(r0 * t);
    dim = t - swindow0 + 1;
    badfs = zeros(dim, 1);
    for i = swindow0:t
        badfs(i - swindow0 + 1) = ADF_IC(y(1:i), adflag, det, IC);
    end
    sadf = max(badfs);
    r2 = swindow0:t;
    rw = r2 - swindow0 + 1;
    bsadfs = zeros(1, dim);
    for v = 1:length(r2)
        swindow = swindow0:r2(v);
        r1 = r2(v) - swindow + 1;
        rwadft = zeros(length(swindow), 1);
        for i = 1:length(swindow)
            rwadft(i) = ADF_IC(y(r1(i):r2(v)), adflag, det, IC);
        end
        bsadfs(v) = max(rwadft);
    end
    if parallel
        parfor v = 1:length(r2)
            swindow = swindow0:r2(v);
            r1 = r2(v) - swindow + 1;
            rwadft = zeros(length(swindow), 1);
            for i = 1:length(swindow)
                rwadft(i) = ADF_IC(y(r1(i):r2(v)), adflag, det, IC);
            end
            bsadfs(v) = max(rwadft);
        end
    end
    gsadf = max(bsadfs);
    result = struct('badfs', badfs, 'bsadfs', bsadfs, 'sadf', sadf, 'gsadf', gsadf);
end