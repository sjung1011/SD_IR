% Clean and Create Data for inflation-relateds regimes: 
% Low and High inflation volatility (S3)
clc
clear all
close all

%% -------------------------------- 
% Deviation-from-target gap 
% ---------------------------------
Tm = readtable('VARDATA.xlsx','Sheet','S1_M');
PCEPI_raw = Tm.PCEPI;

%% Compute annualised month-on-month inflation
pi_pce            = 1200 * diff(log(PCEPI_raw));   % first row lost, monthly to annual (100*12)
Tm.pi_pce = [NaN; pi_pce]; % aling with the date 

%% Deviation-from-target gap
target        = 2.0;                           % Fed’s longer-run PCE goal (annual percent)
gap           = abs(Tm.pi_pce - target);              % |π_t – π*|

%% Binary high-volatility regimes
%  (i) Absolute rule: 
S_abs = double(gap < 1);      % 0 = high-vol; 1 = low-vol, 1pp<pi<3pp

%  (ii) Relative rule: “high vol” if gap exceeds sample 75th percentile
p75   = prctile(gap,75);
S_rel = double(gap < p75);        % 1 = low-vol, 0 = high-vol

Tm.S_abs = S_abs;
Tm.S_rel = S_rel;

%% -------------------------------- 
% Rolling using CPI 
% ---------------------------------
% Setting parameters 
k          = 24;        % trailing-window length in months
percentile = 75;        % threshold percentile for sigma_t

%% compute annualised m/m CPI inflation
CPI  = Tm.CPIAUCSL;                           
pi_cpi   = 1200 * diff(log(CPI));                   % annualised 
Tm.pi_CPI = [NaN; pi_cpi];                          % keep July-62 level row

%% 24-month trailing sigma of inflation
% movstd(x,[k-1 0]) = std over {t-k+1,…,t}
sigma = movstd(pi_cpi,[k-1 0],'omitnan');
Tm.sigma_24 = [NaN; sigma]; % aling with the date

%% --------------------------------------------------------------
% Binary state   (1 = low-vol, 0 = high-vol)   “bad things = 0”
% --------------------------------------------------------------
% Threshold from non-NaN sigma values
thr_75 = prctile(Tm.sigma_24(~isnan(Tm.sigma_24)), percentile);

S_sigma           = NaN(height(Tm),1);
good              = ~isnan(Tm.sigma_24);
S_sigma(good)     = double(Tm.sigma_24(good) < thr_75);   % 1 = low-vol

Tm.S_75thr = S_sigma;

%% Restict data window and Save
startDate = datetime(1962,7,1);
endDate   = datetime(2008,6,1);

keep    = (Tm.observation_date >= startDate) & (Tm.observation_date <= endDate);
Tm      = Tm(keep,:);

outTT = Tm(:,{'observation_date','CPIAUCSL','PCEPI','pi_pce',...
    'pi_CPI','sigma_24','S_abs','S_rel','S_75thr'});

%% Write to Excel
writetable(outTT,'Documents\VARDATA.xlsx','Sheet','S3_Out','WriteMode', 'overwritesheet');

