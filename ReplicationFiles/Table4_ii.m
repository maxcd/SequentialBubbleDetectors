
clear all; clc;
addpath('functions')
%%
m       = 10000; % mc Repititions
rho     = 1.05;
lam0seq = [1:-0.2:0.2];
nlam    = length(lam0seq);
T       = 100;
allow_neg_bubbles = 0;
%% secification of the DGP
mu0     = 1; y00      = 0; beta    = 0; sig2    =1;  r1= 1-lam0seq; DFspec='c'; detrend=0; LMspec='c';
twosided = allow_neg_bubbles;
draws = nan(m,3);

% Boundary b and critival value for CUSUM
lam_cus = 0.919;    
lam_ecus= 2.2;
lam_ocus = 1.902;
r = (1:T-1)'/T;
b = 1+2*r;
rr = (T:-1:1)';


TestNames = {'supDF', 'CUSUM', 'olsCUSUM', 'mCUSUM', 'bCUSUM', 'wCUSUM', 'AHLT$_{10}$'};
ntests = length(TestNames);
RejectionRates = nan(length(lam0seq), ntests);

%%  Start Loop with seed for replicability
lam1seq = lam0seq+r1;
rng('default')
%rng(1234)
for ii=1:length(lam0seq)
    lam0 = lam0seq(ii);
    lam1 =lam1seq(ii);
    t0 = floor(T*lam0);
    t1 = floor(lam1*T);
    

    n0=t0;  % Random Walk perioden
    n1=T-t0;   % explosive perioden

    rhot = [ones(t0-1, 1); repmat(rho, t1-t0, 1); ones(T-t1, 1)];
    H = speye(T) - sparse(2:T, 1:T-1, rhot, T, T);
    
     rej_sdf=0;rej_cus=0; rej_ocus=0; rej_bcus=0; rej_mcus=0; rej_wcus=0; rej_ahlt=0; 
    

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
        % asymptotic critical values
        rej_mcus   = rej_mcus + mres.rej(2);
        rej_cus    = rej_cus  + cres.rej(2);
        % exact critical values
        % rej_mcus   = rej_mcus + (mres.stat > 2.1921);
        %rej_cus    = rej_cus  + (cres.stat  > 0.9196);
        
        
        %  Backward CUSUM
        yR=[cumsum(dy(rr(2:T)))];
        yR=yR/sig/sqrt(T-1);
        ind=abs(yR)>lam_cus*b;
        if sum(ind)>0
           rej_bcus=rej_bcus+1;
        end
        
        
        % OLS CUSUM
        sy=(dy)/sig;
        sy=cumsum(sy)/sqrt(T);
        ostat = abs(max(sy));
        if ostat>lam_ocus
           rej_ocus=rej_ocus+1;
        end
        
        
        %  EoS original (AHLT_10)
        rej_ahlt = rej_ahlt + (EoS_test(y,10,0,twosided) < 0.05);
        
        % weighted-CUSUM
        wres = wCUSUMtest(yges, detrend, twosided);
        rej_wcus=rej_wcus+ wres.rej(2);
        %rej_wcus=rej_wcus+ (wres.stat > 2.12);
        
        % save draws of the CUSUM tests to determine critical values
        draws(i,:) = [cres.stat, mres.stat, wres.stat];
        
        
        % sup DF
        rej_sdf = rej_sdf + supDF(y,detrend);
    
    end
   RejectionRates(ii,:) =  [rej_sdf/m, rej_cus/m ,rej_ocus/m, rej_mcus/m,rej_bcus/m, rej_wcus/m, rej_ahlt/m];
end

%% Size: T100, rho=1
% 0.0507    0.0506    0.0483    0.0501    0.0537    0.0467
RR = array2table([lam0seq', RejectionRates], 'VariableNames', ['r_e', TestNames]);
RR = RR(:,[1, 2, 3, 5, 7, 8])
%% print latex table
tdata = RR{:,:};
cnames = RR.Properties.VariableNames;
[nrows, ncols] = size(tdata);

%ncols = length(cnames);
%nrows = length(lam0seq)
strdgp = ['&\\multicolumn{5}{c}{$T=%2.0f$, $\\rho=%2.2f$, $y_0=%2.0f$, $\\mu_0=%2.0f$, $\\sigma^2=%2.2f$, $\\beta=%2.0f$ }  \\\\ \n '];
strlabel =[ '\\multicolumn{1}{p{2cm}}{\\centering $%s$} ', repmat('& \\multicolumn{1}{p{2cm}}{\\centering $%s$}', 1, ncols-2), ...
    '& \\multicolumn{1}{p{2cm}}{\\centering $%s$}', '\\\\ \\addlinespace \\cmidrule(lr){1-6} \\addlinespace \n'];
str1 = ['%5.1f ', repmat( '& %5.3f ', 1, ncols-1) ,   ' \\\\ \n'];
str2 = '\\multicolumn{1}{p{9cm}}{}   & [%5.2f]   & [%5.2f] &    &     \\\\ \n';


% print the latex table code
fprintf(strdgp, [T, rho, y00, mu0, sig2, beta])
fprintf(strlabel, cnames{:})
for irow=1:nrows
    fprintf(str1, tdata(irow,:))
end
%%



%% 
