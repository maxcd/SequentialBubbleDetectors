function [CUSUM, mCUSUM, wCUSUM] = BubbleCUSUM(y, trend, twosided, est_breakdate, lvl)


%% INPUTS
% y         - Tx1 array of data
% trend     - {0,1} indicator whether to adjust for a trend
% twosided  - do a two-sided test, default=0
% lvl       - scalar indicating the significance level

%% OUTPUTS
% cstat  - CUSUM sequence
% rej   - {0,1} flag with the test decision of the CUSUM test; 1 = reject, 0 = not reject

% mstat  - max(CUSUM) statistic
% mrej   - {0,1} flag with the test decision of the max(CUSUM) test ; 1 = reject, 0 = not reject

%% start the test
    
    if nargin < 2
        trend = 0;
    end

    if nargin<3
        twosided = 0;
    end
    if nargin < 5 || ~ismember(lvl, [10, 5, 1])
       lvl = 5;%[10, 5, 1];
    end

    if nargin < 4 
       est_breakdate=0;
    end

    
    
    % compute some stuff

    T   = length(y);
    dy  = y(2:T)-y(1:T-1);          %  difference




    sig = std(dy);                  % standard dev

    
    % on-sided asymptotic critical values of CUSUM and mCUSUM in each row
    alpha       = [10, 5, 2.5, 1, 0.5];
    CritVals    = [0.74 0.85 0.95 1.06 1.14;...
                    1.64 1.95 2.24 2.57 2.80];
    
    switch trend
        case 1


        % linear bound for the CUSUM test with trend adjustment 
        r=(1:T-2)'/(T-1);
        b_cusum=1+2*r;

        % adjust for a trend recursively
        dy_adj=zeros(T-2,1);        
        for j=2:T-1
            dy_adj(j-1)=sqrt((j-1)/j)*(dy(j)-mean(dy(1:j-1)));
        end
        
        % recursive sigma
        % s2=cumsum(dy_adj.^2)./(1:T-2)';
        % n0=floor(0.3*T);
        % s2(1:n0)=s2(n0)*ones(n0,1);
        
        % compute the cusum test
        y_adj=cumsum(dy_adj);
        cusum_seq=y_adj./sig/sqrt(T-2);
        %cusum_seq_recsigma=y_adj./s2/sqrt(T-1);
        
        n0 = 2;
    case 0
        
        % recursive sigma
        % s2=cumsum((dy-mean(dy)).^2)./(1:T-1)';
        % n0=floor(0.3*T);
        % s2(1:n0)=s2(n0)*ones(n0,1);

        % linear bound for the CUSUM test
        r=(1:T-1)'/T;
        b_cusum=1+2*r;
        y_adj=y(2:T)-y(1);  
        cusum_seq=y_adj./sig/sqrt(T-1);
        %cusum_seq_recsigma=y_adj./s2/sqrt(T-1);

        %cusum_seq= cusum_seq_recsigma;

        n0 = 1;
    end
        

  

     % reject if the indicator exceeds zero anywhere
   
    if twosided % positive and negative bubbles
        cstat = max(abs(cusum_seq)./b_cusum);
        mstat = max(abs(cusum_seq));
        % adjust significance level for twosided tests
        alpha = alpha*2;
    else % one-sided alternative: positive bubbles
        cstat = max(cusum_seq./b_cusum);
%        cstat = max(cstat);
        mstat = max(cusum_seq);
    end
   
    % get critical values to conventional significance levels
    lvls = [10, 5, 1];
    ilvl = find(lvls==lvl);
    crits = CritVals(:, find(ismember(alpha, lvls)));

    crej     = cstat > crits(1,:);
    if sum(crej)==0
        cpval = '>10%';
    else
         [~, tmp] = find(crej, 1, 'last');
        cpval = ['<', num2str(lvls(tmp)), '%'];
    end

    mrej    = mstat > crits(2,:);

   if sum(mrej)==0
        mpval = '>10%';
    else
         [~, tmp] = find(mrej, 1, 'last');
        mpval = ['<', num2str(lvls(tmp)), '%'];
    end
    
    % save results
    CUSUM.type     = 'cusum';
    CUSUM.stat     = cstat;
    CUSUM.rej      = crej;
    CUSUM.alpha    = lvls;
    CUSUM.twosided = twosided;
    CUSUM.pval     = cpval;
    CUSUM.seq      = cusum_seq;
    CUSUM.bound    = b_cusum;
    CUSUM.crit     = crits(1,:);
    [~, first_rej] = max(CUSUM.seq > CUSUM.crit.*CUSUM.bound, [], 1);
    first_rej(~CUSUM.rej) = nan;
    CUSUM.FirstReject = first_rej+n0;



    mCUSUM.type     = 'mcusum';
    mCUSUM.stat     = mstat;
    mCUSUM.rej      = mrej;
    mCUSUM.alpha    = lvls;
    mCUSUM.twosided = twosided;
    mCUSUM.pval     = mpval;
    mCUSUM.seq      = cusum_seq;
    mCUSUM.crit     = crits(2,:);
    [~, first_rej] = max(mCUSUM.seq  > mCUSUM.crit, [], 1);
    first_rej(~mCUSUM.rej) = nan;
    mCUSUM.FirstReject = first_rej+n0;
    
    if est_breakdate
        % CUSUM
        if crej(2)==1
            [bdate_chow,bdate_ML,stats]=getBreakdate(y(1:CUSUM.FirstReject(2)));
            CUSUM.breakdates = [bdate_chow,bdate_ML];
        else
             CUSUM.breakdates = nan(1, 2);
        end
        
        % mCUSUM
        if mrej(2)==1
            [bdate_chow,bdate_ML,stats]=getBreakdate(y(1:mCUSUM.FirstReject(2)));
            mCUSUM.breakdates = [bdate_chow,bdate_ML];
        else
            mCUSUM.breakdates = nan(1, 2);
        end

    end
        
    % weighted-CUSUM
    if nargout==3
        wCUSUM = wCUSUMtest(y, trend, twosided);
        
        if est_breakdate && (wCUSUM.rej(2)==1)
            [bdate_chow,bdate_ML,stats]=getBreakdate(y(1:wCUSUM.FirstReject(2)));
            wCUSUM.breakdates = [bdate_chow,bdate_ML];
        else
            wCUSUM.breakdates = nan(1, 2);
        end

    end
end

function [bdate1,bdate2,stats]=getBreakdate(y)
    %  bdate1:  max(tstat)/Chow
    %  bdate2:  SSR/max(logLikelihood)

    n=size(y,1);
    %y=y-mean(y);
    y=y-y(1);
    dy=y(2:n)-y(1:n-1);
    
    tstat=nan(n,1);
    uu=tstat;
    uu_R = uu;
    n0 = 2; % miminum window length
    for i=n0+1:n-n0

        % OLS in the second part of the sample
        y2  = y(i-1:n-1);
        dy2 = dy(i-1:n-1);
        XX2 = y2'*y2; 
        XY2 = y2'*dy2;
        b2  =  y2\dy2;
        a = XY2; b = XX2;
        
        % Data in the first part of the sample
        y1  = y(1:i-2);
        dy1 = dy(1:i-2);
        b1  = y1\dy1;
        % Unrestricted resiuals in second part of the sample
        u2 = dy2 - y2*b2;

        % Restricted sigma and RSS: Imposing the H0 in the first part of the sample
        u1_R=dy(1:i-2);%-y(1:i-1)*c/d;
        uu_R(i)=u1_R'*u1_R+u2'*u2;
        
        % unrestricted sigma and RSS
        u1=dy1-y1*b1;%c/d;
        uu(i)=u1'*u1+u2'*u2;

        % compute test statistics
        tstat(i)=a/sqrt(b);
    end

    uu_ML = uu_R;

    bdate1 = find(tstat==max(tstat));%+n0;
    bdate2 = find(uu_ML==min(uu_ML));%+n0;
    stats = [tstat, -uu_ML]; 
end
