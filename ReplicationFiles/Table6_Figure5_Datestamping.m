clear; clc;
addpath("functions")
%%
% Produces one Panel of Table 6. Chose the desired number of observations
% below. 
n = 100;            % # observations
% n=200;
% n = 400;

m=10000;            % number of replications
lambda0 = 0.5;       % unit root sample anteil
lambda1 = 1;
n0 = floor(lambda0*n);  % # Random Walk perioden
n1 = n-n0;              % # explosive preioden
rho=1.05;

% DGP
beta=0; mu=1; y0=1; sig2=1; detrend=0;


allow_neg_bubbles = 0;
twosided = allow_neg_bubbles;
%% Asymptotical critial values for SADF and GSADF test from PWY )2015)
T               = [100, 200, 400, 800, 1600]';
r0              = 0.01+1.8./sqrt(T);
sadf_crit95     = [1.3, 1.4, 1.49, 1.53, 1.57]';
gsadf_crit95    = [2, 2.08, 2.2, 2.34, 2.41]';
sadf_asym95     = [1.37, 1.41, 1.49, 1.51, 1.51]';
gsadf_asym95    = [1.89, 2.01, 2.19, 2.2, 2.3]';
CT = table(T, r0, sadf_crit95, gsadf_crit95, sadf_asym95,gsadf_asym95);

det=1;
adflag = 0;
%% Critical values for datestamping #
rng('default')
qe=0.95;
r0 = CT.r0(CT.T==n); % standard choice
swindow0=floor(r0*n);

% parametrical boundary function from Phillips et al. 2011
cv_badf=log(log(swindow0:n)')/100;
cv_sadf = CT.sadf_crit95(CT.T==n)

thisSheet = ['T', num2str(n),'tm',num2str(swindow0)];
%% to use simulated critical values for the SADF, just uncomment the following block

%[~,sheet_names]=xlsfinfo('results/sadf_cvs.xlsx');
% if any(ismember(sheet_names, thisSheet))
%     t = readtable('results/sadf_cvs.xlsx','Sheet',thisSheet);
% else
    % [cv_sadf,cv_badf]=CV_SADF(n,swindow0,0.95);
    % cv_sadf =  [cv_sadf; nan(n-swindow0,1)];
    % t = table(cv_badf, cv_sadf);
    % writetable(t, 'results/sadf_cvs.xlsx', 'Sheet', thisSheet);        
%end
%cv_sadf = t.cv_sadf(1); cv_badf=t.cv_badf;

% only simulate BSADF critical values, if needed (it takes a while)
[~,sheet_names]=xlsfinfo('Critical_Values/gsadf_cvs.xlsx');
if any(ismember(sheet_names, thisSheet))
    t = readtable('Critical_Values//gsadf_cvs.xlsx','Sheet',thisSheet);
else
    [cv_gsadf,cv_bsadf]=CV_GSADF(n,swindow0,0.95);
    cv_gsadf =  [cv_gsadf; nan(n-swindow0,1)];
    t = table(cv_bsadf, cv_gsadf);
    writetable(t, 'Critical_Values//gsadf_cvs.xlsx', 'Sheet', thisSheet);        
end
cv_gsadf = t.cv_gsadf(1); cv_bsadf=t.cv_bsadf;
%%




stats = nan(m,1);
rej_sadf=0; rej_bsadf=0; rej_cus=0; rej_mcus=0; rej_wcus=0; rej6=0;
count=0;
dates_BADF          = nan(m,4);
dates_mC            = nan(m,4);
dates_wC            = nan(m,4);
dates_C             = nan(m,4);
dates_fullSample     = nan(m,4);
rng('default')
for i=1:m

    y = simSingleBubble(n, rho, sig2, lambda0, lambda1, y0, mu, beta, [], allow_neg_bubbles);
    te = floor(n*lambda0)+1;

    % compute a few things
    y_adj=y(2:n)-y(1);          % adjust for starting value
    sig=std(y(2:n)-y(1:n-1));   % Estimate sigma from full sample 
    dy= y(2:n)-y(1:n-1);
    dy_adj=zeros(n-2,1);        % adjust for a trend recursively

    
    % % % GSADF test and BSADF for break date
    res2 = GSADF(y, adflag, det, r0);
    bdate5= NaN;
    %cv_gsadf
    if  res2.gsadf > CT.gsadf_crit95(CT.T==n) %3.3           
        rej_bsadf    = rej_bsadf + 1;
        bdate5   = find(res2.bsadfs > cv_bsadf , 1) + swindow0-1;
    end

    % SDF test and BADF for breakdate
    res1 = SADF(y, adflag, det, r0);
    bdate1 = NaN; bdate1b = NaN;
    if res1.sadf >  CT.sadf_crit95(CT.T==n) %2.5;cv_sadf
        rej_sadf = rej_sadf + 1;
        tmp = res1.badfs(2:end) > cv_badf(2:end);
        bdate1 = find(tmp, 1, "first") + swindow0-1;
        %bdate12 = max(cumsum(~tmp))+1+swindow0-1;
        first_reject = find(res1.badfs > cv_sadf , 1) + swindow0-1;
        %bdate1b = find(res1.badfs > cv_badf , 1) + swindow0;
        n_rej=first_reject;
        [bdate2,bdate3,bdate4]=breakdate_v3b(y(1:n_rej));
        dates_BADF(i,:) = [bdate1, bdate2, bdate3, bdate5];
    end
       

    % CUSUM tests and Chow break date
    [cres, mres, wres] = BubbleCUSUM(y, detrend, twosided, 1);
    % save the date :)
    if cres.rej(2)
        dates_C(i,:) =  [bdate1, cres.breakdates, bdate5];
        rej_cus     = rej_cus  + cres.rej(2);
    end

    if mres.rej(2)
        dates_mC(i,:) = [bdate1, mres.breakdates,  bdate5];
        rej_mcus    = rej_mcus + mres.rej(2);
    end

    if wres.rej(2)
        dates_wC(i,:) = [bdate1, wres.breakdates,  bdate5];
        rej_wcus    = rej_wcus + wres.rej(2);   
    end
    
   
    
end

%

%%
table_dates = dates_wC;

std(table_dates , 'omitmissing')
bine = 1:n/50:n;
bin = discretize(table_dates ,bine);
modebin = mode(bin);
bine = repmat(bine', 1,4 );
%mode_est = 0.5.*[diag(bine(modebin,:))' + diag(bine(modebin-1,:))'];

dateNames = {'SADF',  'Chow', 'ML', 'BSADF'};
tabledates =  dates_wC;
mean_est = round(mean(table_dates , 'omitmissing'), 0);
mode_est = mode(table_dates , 1);
pct_10 = mean((tabledates < n0+1+0.1*n) .* (tabledates> n0+1-0.1*n), 'omitmissing').*100
RMSE= sqrt(mean((tabledates - (n0+1)).^2, 1, 'omitmissing'));
dateStd = std(tabledates, 'omitmissing');
res_tt = array2table([round(mode_est, 0); round(mean_est, 0); round(dateStd,0) ; round(RMSE, 2); round(pct_10, 2)], 'RowNames' , {'Mode', 'Mean', 'Std',  'RMSE',  'P(r^e-0.1<\hat r^e<r^e-0.1)'},...
    'VariableNames', dateNames);

plotlist = [1,4,2];
TableSix = res_tt(:,plotlist)

%% Figure5: Plot Breakdates used with Sample cutting by wCUSUM
plotdates= dates_wC;

textwidth = 426/72.27;
fig1 = figure;
figSize = [textwidth textwidth/4];%
set(fig1, 'PaperUnits', 'inches');
set(fig1, 'Units','inches');
set(fig1, 'PaperSize', figSize);
set(fig1, 'PaperPositionMode', 'auto');
set(fig1, 'Position', [0 0 figSize(1) figSize(2)])
iii = 1;

for iiplot = plotlist
    
    ncols = ceil(sqrt(length(plotlist)))+1;
    nrows = floor(sqrt(length(plotlist)));
    subplot(nrows, ncols, iii)
    hold on;
    
    %binedges = 1:4:n;
    binedges =  (2:n)-0.5;%
    binwidth = 1;
    %binedges = [fliplr(n0+1-2:-4:1), n0+1+2:4:n];
    histogram(plotdates(:,iiplot), binedges, 'EdgeAlpha',0)
    xline(n0+1, 'r-')
    hold off
    xlim([1 n])
    ylim([0 800])
    title(dateNames{iiplot},'Interpreter','latex')
    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
    yaxisproperties= get(gca, 'YAxis');
    yaxisproperties.TickLabelInterpreter = 'latex';   % latex for y-axis
    %if iii>nrows
        xlabel('t', 'Interpreter','latex')
   % end
    box on
    iii=iii+1;
end
tightfig;
