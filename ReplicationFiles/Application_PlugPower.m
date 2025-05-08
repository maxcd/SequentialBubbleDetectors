clear; clc;
addpath("functions")
%% load Bitcoin data
% btc =readtimetable("data\CGS.csv");
% yday = btc.l_AdjClose;
% %figure; plot(btc.obs, yday)
%%
data = readtimetable("data\Plug.csv");
yweek = log(data.Adj_Close);
dates = data.Date;
%figure; plot(btcw.Date, yweek)



%%
y = yweek;
dates = data.Date;
n = length(y);                % # Beobachtungen
% Plug Power. Der Zeitraum ist 2018 - 2021 (Ende Januar ).
sel = logical((dates < datetime('2021-Feb-01')) .* (year(dates) > 2017)) ;
y = y(sel);
dates = dates(sel);
%figure; plot(dates, y)
n=length(y);

%% settings for some tests
% SADF and GSADF
det=1;
adflag = 0;
r0 = 0.01+1.8./sqrt(n); % standard choice
detrend = 0;
twosided = 0;
detrend_adf = 0;
% 0  - Original PWY specification: no detrendinbg but constant in test equation
% 1  - OLS detrending for each ADF window seperately, no constant in test equation
% 2  - recursive detrending, then no constant in test equation.  
%% Critical values for datestamping #
rng('default')
qe=0.95;

swindow0=floor(r0*n);


switch detrend_adf

    case 0 % no detrending, but constant in the Test equation
        thisSheet = ['T', num2str(n),'tm',num2str(swindow0)];
        [~,sheet_names]=xlsfinfo('Critical_Values/sadf_cvs.xlsx');
        if any(ismember(sheet_names, thisSheet))
            t = readtable('Critical_Values/sadf_cvs.xlsx','Sheet',thisSheet);
        else
            [cv_sadf,cv_badf]=CV_SADF(n,swindow0,0.95);
            cv_sadf =  [cv_sadf; nan(n-swindow0,1)];
            t = table(cv_badf, cv_sadf);
            writetable(t, 'Critical_Values/sadf_cvs.xlsx', 'Sheet', thisSheet);        
        end
        

    case 1 % OLS detrending in every Test window
        thisSheet = ['T', num2str(n),'tm',num2str(swindow0)];
        [~,sheet_names]=xlsfinfo('Critical_Values/sadf_demean_cvs.xlsx');
        if any(ismember(sheet_names, thisSheet))
            t = readtable('Critical_Values/sadf_demean_cvs.xlsx','Sheet',thisSheet);
        else
            [cv_sadf,cv_badf]=CV_SDF_demean(n,swindow0,0.95);
            cv_sadf =  [cv_sadf; nan(n-swindow0,1)];
            t = table(cv_badf, cv_sadf);
            writetable(t, 'Critical_Values/sadf_demean_cvs.xlsx', 'Sheet', thisSheet);        
        end

        case 2 % recursive detrending
        thisSheet = ['T', num2str(n),'tm',num2str(swindow0)];
        [~,sheet_names]=xlsfinfo('results/sadf_demean_rec_cvs.xlsx');
        if any(ismember(sheet_names, thisSheet))
            t = readtable('Critical_Values/sadf_demean_rec_cvs.xlsx','Sheet',thisSheet);
        else
            [cv_sadf,cv_badf]=CV_SDF_demean_rec(n,swindow0,0.95);
            cv_sadf =  [cv_sadf; nan(n-swindow0,1)];
            t = table(cv_badf, cv_sadf);
            writetable(t, 'Critical_Values/sadf_demean_rec_cvs.xlsx', 'Sheet', thisSheet);        
        end
end
cv_sadf = t.cv_sadf(1); cv_badf=t.cv_badf;

%%

m=1;
lam=0.82;%0.919;    % für den kritischen Wert von CUSUM
critval1=0.99;

r=(1:n-2)'/n;
b_cusum=1+2*r;
rr=(n:-1:1)';
%rrB=(B-1:-1:1);
%detrend =1;
%r2 = 1;

stats = nan(m,1);
rej_adf=0; rej2=0; rej3=0; rej4=0; rej5=0; rej6=0;
count=0;
dates_BADF = nan(m,4);
rng('default')
i=1;

% compute a few things
y2  = y(2:n)-y(1);          % adjust for starting value
dy  = y(2:n)-y(1:n-1);      % difference
sig = std(dy);              % standard dev
if detrend
    dy_adj=zeros(n-2,1);        % adjust for a trend recursively
    for j=2:n-1
        dy_adj(j-1)=sqrt((j-1)/j)*(dy(j)-mean(dy(1:j-1)));
    end
    %dy_adj = dy(2:end) -mean(dy);
else
    dy_adj = dy; %y(2:n)-y(1);          % adjust for starting value;
end
y_adj=cumsum(dy_adj);
x=[ones(n,1) (1:n)'];
y_bdate=y-x*(x\y);
y_bdate = y_bdate(3:end);

r=(1:length(dy_adj))'/n;
b_cusum=1+2*r;

% 1. sup DF
[reject, dfstat, adf_seq] = supDF(y,detrend);
rej2 = rej2 + (dfstat> 0.79);

% SDF test and BADF for breakdate
% res1 = SADF(y, adflag, det, r0);
% SDF test and BADF for breakdate
switch detrend_adf
    case 0
        res1 = SADF(y, adflag, 1, r0);
    case 1
        res1 = SDF_demean(y, 1, r0);
    case 2
        res1 = SDF_demean_rec(y, r0);
end
sadf_date = NaN; bdate1b = NaN;
if res1.sadf >  cv_sadf
    rej_adf = rej_adf + 1;
    sadf_date = find(res1.badfs > cv_badf , 1) + swindow0-1;
    if detrend_adf==2
       sadf_date=  sadf_date +2;
    end
   
    % find all crossings
    tmp = res1.badfs > cv_badf;
    tmp = tmp(2:end) - tmp(1:end-1);
    sadf_e = find(tmp==1)+ swindow0-1;
    sadf_f = find(tmp==-1)+ swindow0-1;
    sadf_date = sadf_e(end);
end

%% bubble testing with CUSUM variants
% CUSUM, mCUSUM and wCUSUM  tests

% CUSUM tests
[CUSUM, mCUSUM, wCUSUM] = BubbleCUSUM(y, detrend, twosided, 1);
% save the date :)
CUSUM.breakdates
mCUSUM.breakdates
wCUSUM.breakdates


dates_wCUSUM = [sadf_date, wCUSUM.breakdates(1)]
Stats = [res1.sadf, CUSUM.stat, mCUSUM.stat, wCUSUM.stat]
CritVals = [cv_sadf, CUSUM.crit(2), mCUSUM.crit(2), wCUSUM.crit(2)];
%% Store Test Results
testNames = {'supADF','CUSUM', 'mCUSUM', 'wCUSUM'};
Reject= [rej_adf CUSUM.rej(2) mCUSUM.rej(2) wCUSUM.rej(2)]
Reject = array2table([round(Stats, 2); round(CritVals, 2); Reject] ,...
    'VariableNames',  testNames, "RowNames",{'Statistic', 'Critival Value', 'Reject?'})

%% Store Datestamping Results
dateNames = {'PWY', 'Chow'};
BreakDates = {};%:%nan(2, 4);
BreakDates(1,~isnan(dates_wCUSUM')) = cellstr(datestr(dates(dates_wCUSUM(~isnan(dates_wCUSUM')))))'
BreakDates = cell2table(BreakDates, "RowNames",{'Dates'}, "VariableNames",dateNames)


%% Figure 6a: Plot Bubble Statistics

colors = [0.4980    0.2353    0.5529,
    0.0667    0.6471    0.4745,
    0.9490    0.7176    0.0039,
    0.9059    0.2471    0.4549,
    0.5020    0.7294    0.3529,
    0.9020    0.5137    0.0627,
         0    0.5255    0.5843,
    0.8118    0.1098    0.5647,
    0.9765    0.4824    0.4471];

testNames = {'ADF', 'supDF', 'CUSUM', 'mCUSUM', 'wCUSUM'};

    badfs  = res1.badfs;


tdata = [ [nan(floor(r0*n-1), 1); badfs./cv_sadf], ...
    [nan(floor(0.1*n), 1); adf_seq./0.79],...
    [nan(detrend+1, 1); CUSUM.seq./CUSUM.bound./CUSUM.crit(2)],...
    [nan(detrend+1, 1); mCUSUM.seq./mCUSUM.crit(2)],...
    [nan(detrend+1, 1); wCUSUM.seq./wCUSUM.crit(2)]];

plottests =[1,3, 4, 5];

textwidth = 426/72.27;
fig1 = figure;
figSize = [textwidth textwidth/2];%
set(fig1, 'PaperUnits', 'inches');
set(fig1, 'Units','inches');
set(fig1, 'PaperSize', figSize);
set(fig1, 'PaperPositionMode', 'auto');
set(fig1, 'Position', [0 0 figSize(1) figSize(2)])



hold on;

i=1;

colors(3,:) = [];
linespec = {'--', '-', '-.', '--', '-'};
pEoS_stat= {};
for ii=plottests
    pEoS_stat{i} = plot(dates, tdata(:,ii), linespec{i},'color', colors(i,:) ,'LineWidth', 1.2);
    i=i+1;
end

yline(1, 'k:');
hold off

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % latex for y-axis
xlim([dates(1), dates(end)])
%xlabel('t', 'Interpreter','latex')

box on
legend(testNames{plottests}, 'Interpreter','latex', 'location', 'southeast')
legend boxoff



%% Figure 6b: Price Level
textwidth = 426/72.27;
fig1 = figure;
figSize = [textwidth textwidth/2];%
set(fig1, 'PaperUnits', 'inches');
set(fig1, 'Units','inches');
set(fig1, 'PaperSize', figSize);
set(fig1, 'PaperPositionMode', 'auto');
set(fig1, 'Position', [0 0 figSize(1) figSize(2)])

hold on;


linespec = {'--', '-', '-.'};
plot(dates,data.Adj_Close(sel), 'LineWidth',1.2)
offset = [1,-13, -6];
height = [65, 45, 55]
for i=1:2
    if ~isnan(dates_wCUSUM(i))
        b = xline(dates(dates_wCUSUM(i)), linespec{i}, 'color', colors(i,:), 'Linewidth', 1.2);        
        text(dates(dates_wCUSUM(i)+offset(i)), height(i), dateNames{i}, 'color', colors(i,:))
        
    end
end
hold off;

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % latex for y-axis
yaxisproperties.Exponent = 0
yaxisproperties.TickLabelFormat = '%,.0f'
%xlabel('t', 'Interpreter','latex')
ylabel('\$', 'Interpreter', 'latex', 'rotation', 0)
xlim([dates(1), dates(end)])

box on



