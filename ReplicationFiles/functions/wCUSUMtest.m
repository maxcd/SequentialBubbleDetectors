function [res, stat] = wCUSUMtest(y, trend, twosided, cbar, lvl)
% exponentially weighted forward CUSUM test with constant boundary

%% INPUTS
% y         - Tx1 array of data
% trend     - {0,1} indicator whether to adjust for a trend
% twosided  - do a two-sided test, default=0
% cbar      - scalar for computation of the weights using rho=1+cbar/T
% lvl       - default 5% 

%% OUTPUTS
% stat  - eCUSUM statistic
% rej   - {0,1} flag with the test decision; 1 = reject, 0 = not reject 

%% start the test
T = length(y);

if nargin < 2
    trend = 0;
end

if nargin<3
    twosided = 0;
end

if nargin<4
    fixed_weights = 0;
    % choose default value that achieves 50% power in the case that the
    % bubble starts in the middle of the sample
    if trend
        cbar = 2.9;
    else
        cbar = 2;%T*log(rhobar);
    end
elseif cbar == 'fixed'
    fixed_weights = 1;

     if trend
        cbar = 2.9;
    else
        cbar = 2;%2;%T*log(rhobar);
    end
end

if nargin < 5 || ~ismember(lvl, [10, 5, 1])
   lvl = 5;%[10, 5, 1];
end



% compute some stuff
%y   = y(2:T)-y(1);         % eliminate initial value
dy  = y(2:T)-y(1:T-1);      %  1st difference
sigy = std(dy);


rhobar=1+cbar/T;
alpha       = [10, 5, 2.5, 1, 0.5];
switch trend
    case  1 % with trend

        % compute the exponential weights
        switch fixed_weights
            case 0
                w=exp( (0:T-3)'./(T-2)*cbar);
                w=w/sqrt(w'*w);
            case 1
                w=exp( (0:T-3)'./(T-2)*cbar);
                w=w/w(end);
        end
        
        % adjust for a trend recursively
        dy_adj=zeros(T-2,1);        
        for j=2:T-1
            dy_adj(j-1)=sqrt((j-1)/j)*(dy(j)-mean(dy(1:j-1)));
        end

        %dyR=w.*(dy-mean(dy));
        wdy=w.*dy_adj;
        %dyR=dyR(T-1:-1:1);

        sig_dyR = std(wdy);
        %sig_dy =  std(dy);
        switch fixed_weights
            case 0 
                Swy = cumsum(wdy)./sig_dyR./sqrt(T-1);
            case 1
                Swy = cumsum(wdy)./sigy./ sqrt(w'*w)./sqrt(T-1);
        end


        % asymptotical critival values for different cbar
        % first row corresponds to the value of cbar
        CritVals =              [0.5 1.62	1.96	2.25	2.55	2.76;...
                    2   1.63	1.95	2.23	2.59	2.81;...
                    2.9 1.64	1.95	2.23	2.57	2.84;...
                    3   1.64	1.95	2.23	2.57	2.85;...
                    5   1.63	1.96	2.24	2.61	2.83;...
                    10  1.64	1.96	2.24	2.61	2.8];

        
        cseq        = CritVals(:,1);
        CritVals    = CritVals(:,2:end);
        crits       = nan(1, 5);
        for ii=1:5
            if cbar<=cseq(end)
                crits(ii) = interp1(cseq, CritVals(:,ii) ,cbar);
            else
                crits(ii) = CritVals(end);
            end
        end

        n0=2;
    case 0 % no trend

        % compute the exponential weight
        
        r = (0:T-2)'./(T-1);
        w=exp(r*cbar);
        w=w/sqrt(w'*w);
        %v=(exp(2*cbar)-1)/2/cbar;
        %w=w./sqrt(v)/sqrt(T);
        
        
        wdy=w.*dy;
        % sig_dyR = std(wdy);
        % yR=cumsum(wdy)./sig_dyR./sqrt(T-1);
        
        switch fixed_weights
            case 0 
                psi = std(wdy);
                %psi = sqrt(dyR'*dyR);
            case 1
                psi = sigy*w(end);
        end

        %ww2 = (1-1/rhobar^(2*T))/(1-1/rhobar^2);

        Swy = cumsum(wdy)./psi./sqrt(T-1);

        %y_adj=y(2:T)-y(1);  
        %cusum_seq=y_adj/sigy/sqrt(T-1);
        
        % critical values do not depend on cbar without of detrending
        crits = [1.64 1.95 2.24 2.57 2.80];
        n0=1;
end

    if twosided % positive and negative bubbles
        seq = abs(Swy);
        stat = max(abs(Swy));
        % adjust significance level for twosided tests
        alpha = alpha*2;
    else % one-sided alternative: positive bubbles
        seq=Swy;
        stat = max(Swy);
    end
   
    % get critical values to conventional significance levels
    lvls = [10, 5, 1];
    crits = crits(:, find(ismember(alpha, lvls)));

    % determine test decision at conventional levels
    rej     = stat > crits(1,:);
    if sum(rej)==0
        pval = '>10%';
    else
        [~, tmp] = find(rej, 1, 'last');
        pval = ['<', num2str(lvls(tmp)), '%'];
    end

    res.type    = 'wcusum';
    res.cbar    = cbar;
    res.stat    = stat;
    res.seq     = Swy;
    res.rej     = rej;
    res.crit    = crits;
    res.alpha   = lvls;
    res.twosided= twosided;
    res.pval    = pval;
    res.seq     = seq;
    [~, first_rej] = max(res.seq  > res.crit, [], 1);
    first_rej(~res.rej) = nan;
    res.FirstReject = first_rej+n0;
end