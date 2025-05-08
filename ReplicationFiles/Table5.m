clear; clc;
addpath('functions')

% Set simulation Values from PSY
%%
DGP     = 0;
% 0  - constant mean
% 22 - trend
m       = 1000; % mc Repititions
rho     = 1.05;

%% chose volatility types
Types = {'constant', 'garch', 'ST-up','cos', };%
HetNames = {'Homoscedasticity', 'GARCH(1,1)', 'ST-up', 'cosine'};
% constant
% garch
% cos
% cos-down
% cos-up
% SM-down
% SM-up
lam0seq = [1:-0.2:0.6];
nlam    = length(lam0seq);
T       = 100;

% set seed for replicability
twosided=1;


%% secification of the DGP
switch DGP
    case 0 % DGP 1: no constant, no initial value
        mu0     = 0; y00      = 0; beta    = 0; sig2    = 1;  r1= 1-lam0seq; DFspec='n'; detrend=0; LMspec='n';
        SelectTests = [1, 2, 6, 5, 3];
    case 22
        mu0     = 1; y00      = 0; beta    = 1; sig2    = 1;  r1= 1-lam0seq; DFspec='ct'; detrend=1;
        SelectTests = [1, 2, 6, 5, 4];
end
%%
rng('default')
RR_list = {};
lam1seq = min(lam0seq+r1, 1);
TestNames = {'supDF', 'CUSUM',  'AHLT$_{10}$',  'pEoS$_{0.3}$', 'wCUSUM','mCUSUM'};
ntests = length(TestNames);

TableFive =[];

for ii=1:length(lam0seq)
        
    RejectionRates = nan(length(Types), ntests);

     for itype =1:length(Types)
        VolaType = Types{itype};


        lam0 = lam0seq(ii);
        lam1 =lam1seq(ii);
        t0 = floor(T*lam0);
        t1 = floor(lam1*T);
        
    
        n0=t0;  % Random Walk perioden
        n1=T-t0;   % explosive perioden
    
        rhot = [ones(t0-1, 1); repmat(rho, t1-t0, 1); ones(T-t1, 1)];
        H = speye(T) - sparse(2:T, 1:T-1, rhot, T, T);
        
        rej_sdf=0; rej_cus=0; rej3=0; rej4=0; rej_ahlt=0; 
        rej61=0; rej_eos=0; rej8=0; rej_wcus=0; rej9=0;
        rej_mcus=0;
    
        for i=1:m
            eps = randn(T,1);
    
            [sig2, u, ~] = generateSig2(VolaType, T, eps);
            %u  = sqrt(sig2).*eps;
            det = ones(T,1)*mu0 + (1:T)'*beta;
            y = det+H \([y00; zeros(T-1,1)] + u);
            
            % compute some stuff
            %y   = y(2:T)-y(1);          % eliminate initial value
            dy  = y(2:T)-y(1:T-1);      %  difference
            sig = std(dy);              % standard dev
            dy_adj=zeros(T-2,1);        % adjust for a trend recursively
            for j=2:T-1
                dy_adj(j-1)=sqrt((j-1)/j)*(dy(j)-mean(dy(1:j-1)));
            end
    


            % forward CUSUM and max CUSUM
            [cres, mres] = BubbleCUSUM(y, detrend, twosided);
            % asymptotic critical values
            rej_mcus   = rej_mcus + mres.rej(2);
            rej_cus    = rej_cus  + cres.rej(2);

                   
         
            % eCUSUM
            wres = wCUSUMtest(y, detrend, twosided);
            rej_wcus=rej_wcus+ wres.rej(2);
           
            
            % sup DF
            %rej6 = rej6 + supDF(y,detrend);
            [reject, dfstat, pval] = supDF(y,detrend);
            rej_sdf = rej_sdf + reject;%(dfstat> 0.75);
            
            % AHLT (10) test
            rej_ahlt = rej_ahlt + (EoS_test(y,10,1,1) < 0.05);
    
    
             % EoS(30) parameteric
            k03=floor(0.3*T);
            w0=1.0.^(1:k03)';
            weos=[ -sum(w0)/(T-k03-1)*ones(T-k03-1,1) ; w0 ];
            a=weos.*(dy-mean(dy));
            a=sum(a)/sqrt(a'*a);
            if abs(a)>1.96
                rej_eos=rej_eos+1;
            end
            
        
      
            
         
        
        end
        	

       RejectionRates(itype,:) =  [rej_sdf/m, rej_cus/m ,rej_ahlt/m, rej_eos/m, rej_wcus/m, rej_mcus/m];
       


    end
    
    %% Size: T100, rho=1
    % 0.0507    0.0506    0.0483    0.0501    0.0537    0.0467
    tmp = array2table(HetNames', 'VariableNames', {'vola type'} );
    tmp2 = array2table(ones(length(HetNames'),1)*lam0, 'VariableNames', {'r_e'} );
    RR = [tmp2, tmp, array2table(RejectionRates(:,SelectTests), 'VariableNames', TestNames(SelectTests))];
    TableFive = [TableFive; RR];
    


end

TableFive
