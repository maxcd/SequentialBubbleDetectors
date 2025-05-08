function [CUSUM, mCUSUM, wCUSUM] = BubbleCUSUMmonitoring(y, T0, Tmax, cbar, trend, twosided)
% copute the statistic for the monitoring procedures based on the CUSUM,
% mCUSUM, wCUSUM

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
    
    if nargin < 3

        disp('Monitoring requires at least the first three inputs, includig T0 and Tmax')

    elseif length(y)>(T0+Tmax)
    
        disp('Time Series is longer thatn the Monitoring Horizon!')
    end
    
    if nargin < 5
        trend = 0;
    end

    if nargin<6
        twosided = 0;
    end


    
    % weight for the wCUSUM test
    if nargin<4
        if trend
            cbar = 2.9;
        else
            cbar = 2.1;
        end
    end
    
    

    % compute some stuff
    T   = length(y);
    t   = T-T0;
    dy  = y(2:T)-y(1:T-1);          %  difference
    %sig = std(dy);                  % standard dev


    dy0 = dy(1:T0-1);
    dy = dy(T0:T0+min(Tmax, T-T0)-1);



    sig = std(dy0);
    
    
    switch trend
        case 1 % trend case


        % linear bound for the CUSUM test with trend adjustment 
        r=(1:t-2)'/(Tmax-1);
        b_cusum=1+2*r;

        % adjust for a trend recursively
        dy_adj=zeros(t-2,1);        
        for j=2:t-1
            dy_adj(j-1)=sqrt((j-1)/j)*(dy(j)-mean(dy(1:j-1)));
        end
        
        % compute the cusum test
        y_adj=cumsum(dy_adj);
        cusum_seq=y_adj/sig/sqrt(Tmax-2);



         % weighted CUSUM sequence
        r = (0:t-2)'./(Tmax-2);
        % weights
        w=exp(r*cbar);
        % normalize
        w=w/sqrt(w'*w);
       
        
        % weighted differences
        wdy=w.*dy_adj;
        % under homoscedasticity psi equals sig
        psi = sig;
        Swy = cumsum(wdy)./psi./sqrt(Tmax-1);
        % divide by the window specific bound wt for monitoring
        Swy = Swy./w(end);
        

        case 0 % no trend case
        
        % linear bound for the CUSUM test
        r=(1:t)'/Tmax;
        b_cusum=1+2*r;
        y_adj=cumsum(dy); 
        
        % CUSUM sequence
        cusum_seq=y_adj/sig/sqrt(Tmax-1);

        % weighted CUSUM sequence
        r = (0:t-1)'./(Tmax-1);
        % weights
        w=exp(r*cbar);
        % normalize
        w=w/sqrt(w'*w);
        
        %v=(exp(2*cbar)-1)/2/cbar;
        %w=w./sqrt(v)/sqrt(T);
        
        % weighted differences
        wdy=w.*dy;
        % under homoscedasticity psi equals sig
        psi = sig;
        Swy = cumsum(wdy)./psi./sqrt(Tmax-1);
        % divide by the window specific bound wt for monitoring
        Swy = Swy./w(end);
        


    end
        
    % apply absolute value for a twosided test
    if twosided % positive and negative bubbles
        cstat = max(abs(cusum_seq)./b_cusum);
        mstat = max(abs(cusum_seq));
        Swy = abs(Swy);
        wstat = max(abs(Swy));
        % adjust significance level for twosided tests
        %alpha = alpha*2;
    else % one-sided alternative: positive bubbles
        cstat = max(cusum_seq./b_cusum);
        mstat = max(cusum_seq);
        wstat = max(Swy);
    end
   
   
    % save results
    CUSUM.type     = 'cusum';
    CUSUM.stat     = cstat;
    CUSUM.twosided = twosided;
    CUSUM.seq      = cusum_seq;
    CUSUM.bound    = b_cusum;


    mCUSUM.type     = 'mcusum';
    mCUSUM.stat     = mstat;
    mCUSUM.twosided = twosided;
    mCUSUM.seq      = cusum_seq;


    wCUSUM.type    = 'wcusum';
    wCUSUM.cbar    = cbar;
    wCUSUM.stat    = wstat;
    wCUSUM.seq     = Swy;
    wCUSUM.twosided= twosided;

end

