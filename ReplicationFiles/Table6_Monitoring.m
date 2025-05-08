clear all; clc;
addpath('functions')
%%
% Set simulation Values from PSY
DGPseq     = [2, 3];
% 2 - homoskedastic

PositiveBubblesOnly = 1;
m       = 10000; % mc Repititions

%% DGP specification
rho     = 1.05;
cbar    = 2.1;
lam0seq = [1:-0.2:0.2];%[1:-0.1:0.6];
nlam    = length(lam0seq);
T       = 50;
T0      = 50;
TT      = T+T0;%250;%250;%100
tm      = TT-T;         % Start of monitorin process

% set seed for replicability
rng('default')

het=0;
%% secification of the DGP
RejectionRates = [];
LamEst = [];

% for sDF
swindow0=T0+1;
r0 = swindow0/TT;%CT.r0(CT.T==n); % standard choice

%% DGP
mu0     = 1; y00      = 0; beta    = 0;   r1= 1-lam0seq;  detrend=0; adfdet=1; sig2=1;

%% simulate critical values for the sadf test 
[cv_sadf,cv_badf]=CV_SADF(TT,swindow0,0.95);
%[cv_sadf,cv_badf]=CV_SDF_demean_rec(TT,swindow0,0.95,1);
cv_sadf =  [cv_sadf; nan(TT-swindow0,1)];
t = table(cv_badf, cv_sadf);
cv_sadf = t.cv_sadf(1); cv_badf=t.cv_badf;
   
    
if PositiveBubblesOnly
    allow_neg_bubbles = 0;
    twosided=0;
else
    twosided=1;
    allow_neg_bubbles = 1;
end
    

    
    %[cv_sadf,cv_badf]=CV_SADF(T,swindow0,0.95);
    % parametrical boundary function from Phillips et al. 2011
    %cv_badf_par=log(log(swindow0:T)')/100;
    
% critical values of the wCUSUM test
%Table 2: One-sided critical values for monitoring
%columns cbar 10% 5% 2.5% 1% 0.5%
wCusumCrit = [ 0 1.62 1.95 2.25 2.62 2.88;
1 1.28 1.5 1.69 1.93 2.08;
2 1.10 1.26 1.41 1.59 1.71;
2.1 1.09 1.25 1.39 1.57 1.69;
3 0.98 1.12 1.24 1.38 1.48;
4 0.90 1.01 1.12 1.24 1.33;
5 0.84 0.94 1.03 1.14 1.22;
6 0.79 0.88 0.96 1.06 1.14;
8 0.72 0.80 0.87 0.95 1.0];

idx = find(wCusumCrit(:,1) == cbar);
wCUSUM05 = wCusumCrit(idx,3);
    
    
    %%
    lam1seq = min(lam0seq+r1, 1);
    
    TestNames = {'CUSUM', 'mCUSUM', 'wCUSUM', 'sup DF'};
    DateNames = {'Chow', 'SADF'};
    ModeNames = {'Chow', 'ML', 'SADF (1st)', 'SADF (2nd)'};
    ntests = length(TestNames);
    ndates = length(DateNames);
    RejectionRates_tmp = nan(length(lam0seq), ntests);
    lam0hat_tmp = nan(length(lam0seq), ntests);
    DelayMean = nan(length(lam0seq), ntests);
    RejectTooEarly = nan(length(lam0seq), ntests);

    
    
   for ii=1:length(lam0seq)
        
        lam1 =lam1seq(ii);
        te = floor(T*lam0seq(ii))+T0;
        t1 = floor(lam1*T)+T0;
        lam0 = te/TT; % for simulation of the DGP
    
        n0=te;      % Random Walk perioden
        n1=TT-te;   % explosive perioden
    
        rhot = [ones(te-1, 1); repmat(rho, t1-te, 1); ones(TT-t1, 1)];
        H = speye(TT) - sparse(2:TT, 1:TT-1, rhot, TT, TT);
        
        rej2=0; rej1=0; rej3=0; rej4=0;
        rej_cusum = zeros(1, 3);
        n_rej2=0; n_rej1=0; n_rej3=0; n_rej4=0;
        
        store_dates     = nan(m, ndates);
        store_detection = nan(m, 4);
        store_delay     = nan(m, 4);
        for i=1:m
            n_rej2=nan; n_rej1=nan; n_rej3=nan; n_rej4=nan;
            date_chow = nan; date_ml=nan;
            % simulate data
            
            y = simSingleBubble(TT, rho, sig2, lam0, lam1, y00, mu0, beta, [], allow_neg_bubbles);
            

            % CUSUM
            stats_seq   = nan(TT, 3);
            cus_end     = nan(TT, 1);
            %rej_seq   = nan(TT, 3);
            % estimate sigma
            dy = y(2:end) - y(1:end-1);
            sighat = std(dy(1:T0-1));
            
            for it=T0+2:TT
                [cres, mres, wres] = BubbleCUSUMmonitoring(y(1:it), T0, T, cbar);
                %[cres, mres, wres] = BubbleCUSUMmonitoring(y(1:it), detrend,  twosided, cbar, sighat);
                % asymptotic critical values
                 stats_seq(it,:) = [cres.stat, mres.stat, wres.stat];
                 cus_end(it) = cres.seq(end);
                 %rej_seq(it,:) = [crcres.states.rej(2), mres.rej(2), wres.rej(2)];
            end
            



           Reject = stats_seq >  [0.85   1.96    wCUSUM05];


            any_rej = any(Reject);
            rej_cusum = rej_cusum + any_rej;
            
            % date-stamping (not reported in the paper)
            FirstReject = nan(1,3);
            if any_rej(3) % if the wCUSUM rejects anywhere

                [~, FirstReject] = max(Reject , [], 1);
                FirstReject(~any_rej) = nan;
                dy=y(2:TT)-y(1:TT-1);
                 n2=FirstReject(3);
                 switch detrend
                     case 0
                        ybdate = y;
                        loose_obs = 0;
                     case 1
                        x=[ones(TT,1) (1:TT)'];
                        ybdate=y-x*(x\y);
                        dy_adj=zeros(TT-2,1);        % adjust for a trend recursively
                        for j=2:TT-1
                            dy_adj(j-1)=sqrt((j-1)/j)*(dy(j)-mean(dy(1:j-1)));
                        end
                        ybdate=[zeros(2, 1) ;cumsum(dy_adj)];
                        loose_obs = 0;                   
                 end
                
                % estimate the break date
                [date_chow, ChowStat]=getChowBreakdate(ybdate(1:n2-loose_obs), T0);
                
            end

            
            % sup DF
            res4 = SADF(y, 0, 1, r0);

            date_sadf = NaN;
            detect_sadf = NaN;
            if res4.sadf >  cv_sadf %CT.sadf_crit95(CT.T==n) %2.5;cv_sadf
                rej4 = rej4 + 1;
                % date-stamping
                date_sadf = find(res4.badfs > cv_badf , 1) + swindow0-1;
                % detection delay
                detect_sadf = find(res4.badfs > cv_sadf , 1) + swindow0-1;
            end
            dates_tmp = [date_chow, date_sadf];
            store_dates(i,:) = dates_tmp;
            store_delay(i,:) = [FirstReject, detect_sadf]-(te+1);
            store_detection(i,:) = [FirstReject, detect_sadf];

        end

       % compute the delay 
       posDelay = store_delay>0;
       RejectTooEarly(ii,:) = mean(~posDelay);
       store_delay(~posDelay) = 0;

       % results for the table
       DelayMean(ii,:)      = mean(store_delay, 1, 'omitmissing');
       RejectionRates_tmp(ii,:) =  [rej_cusum/m rej4/m];
  end
    

Rejections = array2table(RejectionRates_tmp, 'VariableNames', TestNames);
DelayMean = array2table(DelayMean, 'VariableNames', TestNames);
     
            
    
%% Compute Results

% date stamping
lamtab = array2table(lam0seq', 'VariableNames', {'r_e'});
breaktab = array2table(min(floor(lam0seq'.*T)+1,TT), 'VariableNames', {'r_e'});
TableSix = table(breaktab, round(Rejections, 3), round(DelayMean, 0))

    



