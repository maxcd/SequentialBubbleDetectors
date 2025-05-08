clear all; clc;
addpath('functions')



m       = 10000; % mc Repititions
rho     = 1.05;
lam0seq = [1:-0.2:0.2];
nlam    = length(lam0seq);
T       = 100;



for subtabs=[1:2]

    switch subtabs
        case 1 
            tabtitle = 'Table3_i_twosided',
            allow_neg_bubbles = 1;
        case 2 
            tabtitle = 'Table3_ii_onesided';
            allow_neg_bubbles = 0;
    end




%% secification of the DGP
mu0     = 1; y00      = 0; beta    = 0; sig2    =1;  r1= 1-lam0seq;  

% specify the tests
detrend=0; 
twosided = allow_neg_bubbles;

% storage matrix
TestNames   = {'supDF', 'CUSUM', 'mCUSUM', 'wCUSUM', 'AHLT$_{10}$'};
ntests      = length(TestNames);
RejectionRates = nan(length(lam0seq), ntests);

%%  Start Loop with seed for replicability
lam1seq = lam0seq+r1;
rng('default')

for ii=1:length(lam0seq)
    lam0 = lam0seq(ii);
    lam1 =lam1seq(ii);
    t0 = floor(T*lam0);
    t1 = floor(lam1*T);
   
    
     rej_sdf=0;rej_cus=0; rej_mcus=0; rej_wcus=0; rej_ahlt=0; 
    

    for i=1:m

        % simulate Bubble
        yges = simSingleBubble(T, rho, sig2, lam0, lam1, y00, mu0, beta, [], allow_neg_bubbles);
        
        % de-mean, differencing, sigma
        y=yges(2:T)-yges(1);
        dy = yges(2:T)-yges(1:T-1);
        sig = std(dy);
        
        %  CUSUM
        % forward CUSUM and max CUSUM
        [cres, mres] = BubbleCUSUM(yges, detrend, twosided, 5);
        % store test decision
        rej_mcus   = rej_mcus + mres.rej(2);
        rej_cus    = rej_cus  + cres.rej(2);
        
        %  EoS original (AHLT_10)
        rej_ahlt = rej_ahlt + (EoS_test(y,10,detrend,twosided) < 0.05);
        
        % weighted-CUSUM
        wres = wCUSUMtest(yges, detrend, twosided);
        rej_wcus=rej_wcus+ wres.rej(2);
                
        % sup DF
        rej_sdf = rej_sdf + supDF(y,detrend);
    
    end
   RejectionRates(ii,:) =  [rej_sdf/m, rej_cus/m, rej_mcus/m, rej_wcus/m, rej_ahlt/m];
end

 
RejectionRates = round(RejectionRates, 3);
switch subtabs
    case 1 
        tab3_i_twosided = array2table([lam0seq', RejectionRates], 'VariableNames', ['r_e', TestNames])
    case 2
        tab3_ii_onesided = array2table([lam0seq', RejectionRates], 'VariableNames', ['r_e', TestNames])
end

end