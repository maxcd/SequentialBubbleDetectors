function [tehat, ChowStat]=getChowBreakdate(y, T0)
    %  tehat: bubble emeergence estimation with  max(tstat)/max(Chow)

    n=size(y,1);

    y=y-y(1);
    dy=y(2:n)-y(1:n-1);
    
    tstat=nan(n,1);

    nmin     = 2; % miminum window length

    if nargin<2
        T0 = nmin;
    end

    for i=T0+1:n-nmin

        % OLS in the second part of the sample
        y2  = y(i-1:n-1);
        dy2 = dy(i-1:n-1);
        XX2 = y2'*y2; 
        XY2 = y2'*dy2;
        a = XY2; b = XX2;
        

        % compute test statistics
        tstat(i)=a/sqrt(b);
    end


    tehat = find(tstat==max(tstat));%+n0;
    ChowStat = tstat; 
end