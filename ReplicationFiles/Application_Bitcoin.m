clear all; clc;
addpath("functions")
%% load Bitcoin data
data = readtimetable("data\BTC_USD_Dec24.csv");
data = retime(data, 'weekly', 'lastvalue');
yweek = log(data.adj_close);
%figure; plot(data.Date, yweek)
%%
y = yweek;
dates = data.date;
n = length(y);                % # Beobachtungen

sel = logical((year(dates) < 2025).* (dates > datetime('2022-09-30', 'InputFormat', 'uuuu-MM-dd'))) ; %2021 % dates> '2023-Jul-01'
%sel = logical((year(dates) < 2025).* (dates>['2021-Dec-31']));
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
detrend_adf = 0; 
twosided = 0;
% 0  - Original PWY specification: no detrendinbg but constant in test equation
% 1  - OLS detrending
% 2  - recursive detrending, then no constant in test equation.  
%% Critical values for datestamping #
rng('default')
qe=0.95;

swindow0=floor(r0*n);

switch detrend_adf

    case 2 % recursive detrending
        thisSheet = ['T', num2str(n),'tm',num2str(swindow0)];
        [~,sheet_names]=xlsfinfo('Critical_Values/sadf_demean_rec_cvs.xlsx');
        if any(ismember(sheet_names, thisSheet))
            t = readtable('Critical_Valuess/sadf_demean_rec_cvs.xlsx','Sheet',thisSheet);
        else
            [cv_sadf,cv_badf]=CV_SDF_demean_rec(n,swindow0,0.95);
            cv_sadf =  [cv_sadf; nan(n-swindow0,1)];
            t = table(cv_badf, cv_sadf);
            writetable(t, 'Critical_Values/sadf_demean_rec_cvs.xlsx', 'Sheet', thisSheet);        
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
end
cv_sadf = t.cv_sadf(1); cv_badf=t.cv_badf;

%%
rej_adf=0; 

rng('default')
i=1;

% compute a few things
y2  = y(2:n)-y(1);          % adjust for starting value
dy  = y(2:n)-y(1:n-1);      % difference

if detrend
    dy_adj=zeros(n-2,1);        % adjust for a trend recursively
    for j=2:n-1
        dy_adj(j-1)=sqrt((j-1)/j)*(dy(j)-mean(dy(1:j-1)));
    end
else
    dy_adj = dy; %y(2:n)-y(1);          % adjust for starting value;
end


%% SDF test and BADF for breakdate
switch detrend_adf
    case 1
        res1 = SDF_demean(y, 1, r0);
    case 0
        res1 = SADF(y, adflag, 1, r0);
    case 2
        res1 = SDF_demean_rec(y, r0);
end
sadf_date = NaN; bdate1b = NaN;
if res1.sadf >  cv_sadf
    rej_adf = rej_adf + 1;
    sadf_date = find(res1.badfs(2:end) > cv_badf(2:end) , 1) + swindow0-1;
    if detrend_adf==2
       sadf_date=  sadf_date +2;
    end
    %bdate1b = find(res1.badfs > cv_badf , 1) + swindow0;
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




%% date-stamping summary
est_dates = [sadf_date, wCUSUM.breakdates(1)]
Stats = [res1.sadf, CUSUM.stat, mCUSUM.stat, wCUSUM.stat];
CritVals = [cv_sadf, CUSUM.crit(2), mCUSUM.crit(2), wCUSUM.crit(2)];
%% Store Test Results
testNames = {'supADF','CUSUM', 'mCUSUM', 'wCUSUM'};
Reject= [rej_adf CUSUM.rej(2) mCUSUM.rej(2) wCUSUM.rej(2)]
Reject = array2table([round(Stats, 2); round(CritVals, 2); Reject] ,...
    'VariableNames',  testNames, "RowNames",{'Statistic', 'Critival Value', 'Reject?'})

%% Store Datestamping Results
dateNames = {'PWY', 'Chow'};
BreakDates = {};%:%nan(2, 4);
BreakDates(1,~isnan(est_dates')) = cellstr(datestr(dates(est_dates(~isnan(est_dates')))))'
BreakDates = cell2table(BreakDates, "RowNames",{'Dates'}, "VariableNames",dateNames)
%% Figure 7b): Normalized bubble statistics

testNames = {'ADF', 'CUSUM', 'mCUSUM', 'wCUSUM'};

  badfs  = res1.badfs;
tdata = [ [nan(floor(r0*n-1), 1); badfs./cv_sadf], ...
    [nan(detrend+1, 1); CUSUM.seq./CUSUM.bound./CUSUM.crit(2)],...
    [nan(detrend+1, 1); mCUSUM.seq./mCUSUM.crit(2)],...
    [nan(detrend+1, 1); wCUSUM.seq./wCUSUM.crit(2)]];

plottests =[1, 2, 3, 4];

textwidth = 426/72.27;
fig1 = figure;
figSize = [textwidth textwidth/2];%
set(fig1, 'PaperUnits', 'inches');
set(fig1, 'Units','inches');
set(fig1, 'PaperSize', figSize);
set(fig1, 'PaperPositionMode', 'auto');
set(fig1, 'Position', [0 0 figSize(1) figSize(2)])


%ylim([-4, 2])
hold on;

i=1;
colors = [0.4980    0.2353    0.5529,
    0.0667    0.6471    0.4745,
    0.9490    0.7176    0.0039,
    0.9059    0.2471    0.4549,
    0.5020    0.7294    0.3529,
    0.9020    0.5137    0.0627,
         0    0.5255    0.5843,
    0.8118    0.1098    0.5647,
    0.9765    0.4824    0.4471];
colors(3,:) = [];
linespec = {'--', '-', '-.', '--'};
a= {};
for ii=plottests
    a{i} = plot(dates, tdata(:,ii), linespec{i},'color', colors(i,:) ,'LineWidth', 1.2);
    i=i+1;
end

yline(1, 'k:');


hold off

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % latex for y-axis
xlim([dates(1), dates(end)])
ylim([-1.7 1.5])
%xlabel('t', 'Interpreter','latex')

box on
legend(testNames{plottests}, 'Interpreter','latex', 'location', 'southeast')
legend boxoff


%% Fibure 7a: Data and estimated bubble emergence
textwidth = 426/72.27;
fig1 = figure;
figSize = [textwidth textwidth/2];%
set(fig1, 'PaperUnits', 'inches');
set(fig1, 'Units','inches');
set(fig1, 'PaperSize', figSize);
set(fig1, 'PaperPositionMode', 'auto');
set(fig1, 'Position', [0 0 figSize(1) figSize(2)])

hold on;

plot(dates,data.adj_close(sel), 'LineWidth',1.2);


offset = [0,-10,-2];
for i=1:2
    if ~isnan(est_dates(i))
        b = xline(dates(est_dates(i)), linespec{i}, 'color', colors(i,:),  'LineWidth',1.2);
        text(dates(est_dates(i)+offset(i)), 10000, dateNames{i}, 'color', colors(i,:), 'LineWidth',1.2);
    end
end
hold off;
hold off;

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % latex for y-axis
yaxisproperties.Exponent = 0;
yaxisproperties.TickLabelFormat = '%,.0f';
%xlabel('t', 'Interpreter','latex')
ylabel('\$', 'Interpreter', 'latex', 'rotation', 0)
ylim([0, 110000]);
xlim([dates(1), dates(end)]);
box on
