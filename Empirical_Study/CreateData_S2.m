% Clean and Create Data for inflation-relateds regimes: 
% Anchored and unanchored (S2)

clc
clear all
close all

%% -------------------------------- 
% Michigan Survey (HH) 
% ---------------------------------
% Import the Michigan survey data
Tm = readtable('VARDATA.xlsx','Sheet','S2_M');  
% rename and date vector                                    
iqrHH  = Tm.('InterquartileRange_75th_25th');

if isnumeric(Tm.Month)                        % case 1: Month is 1–12
    monthNum = Tm.Month;
else                                          % case 2: Month is text
    monthNum = month(datetime(Tm.Month,'InputFormat','MMM'));  % or 'MMMM'
end

dates  = datetime(Tm.Year, monthNum, 1);       % first day of each month
Tm.Date = dates; 

%% Define the two sample windows
lastDate    = datetime(2008,6,1);          % hard cutoff at 2008-07-01
cutMask     = dates <= lastDate;           % keep only rows up to July 2008
Tm          = Tm(cutMask,:);                % drop later rows entirely

dStart78   = datetime(1978,1,1);
dStart81   = datetime(1981,7,1);
dEnd       = datetime(2008,6,1);        

mask78 =  dates>=dStart78 & dates<=dEnd;   % 1978-01 – 2008-06
mask81 =  dates>=dStart81 & dates<=dEnd;   % 1981-07 – 2008-06
%% Compute the 40-th percentile of the (non-missing) IQR time-series
data78   = iqrHH(mask78);
data81   = iqrHH(mask81);

p40_78   = prctile(data78(~isnan(data78)),40);   % 1978–2008 window
p40_81   = prctile(data81(~isnan(data81)),40);   % 1981–2008 window

%% Generate the binary regime indicator
S78 = NaN(height(Tm),1);
S81 = NaN(height(Tm),1);

S78(mask78) = double(iqrHH(mask78) <= p40_78);
S81(mask81) = double(iqrHH(mask81) <= p40_81);

%% Attach the indicator to the table 
Tm.S_anchor78 = S78;
Tm.S_anchor81 = S81;

% Quick sanity check (fraction of achored)
% returns a share in the 0.55–0.65 range and that the anchored months cluster in the post-Volcker period. 
% Adjust the percentile (e.g., 35 % or 45 %) only if the split looks implausible.
mean(S78==1,'omitnan') % 0.41
mean(S81==1,'omitnan') % 0.37

%% -------------------------------- 
% SPF (Professionals) one-year-ahead CPI IQR
% ---------------------------------
% Import the SPF data
Tq = readtable('VARDATA.xlsx','Sheet','S2_Q');

% rename and date vector
iqrPRO_Q = Tq.CPI_D1_T_4_; 
datesQ = datetime(Tq.Survey_Date_T_, ...
                  'InputFormat','yyyyQQQ','Format','yyyy-MM-dd');
%% Replicate each quaterly value across the 3 months of it
yr        = year(datesQ);
qr        = quarter(datesQ);
mStart    = (qr-1)*3 + 1;                        % 1, 4, 7, 10
yrsRep     = repelem(yr, 3);                 %  N×1  →  (3N)×1
monthsRep  = repelem(mStart, 3) + ...        %  (3N)×1
             repmat( (0:2)', numel(mStart), 1 );

datesM   = datetime(yrsRep, monthsRep, 1); 
datesM    = datesM(:);                           % column vector: Jul,Aug,Sep,…

iqrPRO_M  = repelem(iqrPRO_Q,3);                 % same length as datesM

%% Data trimming
maskPRO   = datesM >= dStart81 & datesM <= dEnd;
datesM    = datesM(maskPRO);
iqrPRO_M  = iqrPRO_M(maskPRO);
%% 40-th-percentile threshold and binary regime S_PRO
p40_PRO   = prctile(iqrPRO_M(~isnan(iqrPRO_M)),40);
S_PRO     = double(iqrPRO_M <= p40_PRO);         % 1 = anchored
S_PRO(isnan(iqrPRO_M)) = NaN;

%% Align S_PRO to the full Michigan calendar (fill NaNs before 1981-07)
allDates      = Tm.Date;                      % 1978-01 … 2008-07  (N×1)
S_PRO_full    = NaN(size(allDates));          % pre-allocate with NaNs

[isHit,loc]   = ismember(datesM, allDates);   % match SPF months to Michigan dates
S_PRO_full(loc(isHit)) = S_PRO(isHit);        % fill the matched rows

% Quick sanity check (fraction of achored)
mean(S_PRO_full==1,'omitnan') % 0.36
%% Save 
Tm.S_anchorPRO = S_PRO_full;
%% Out 
outTT = Tm(:,{'Date','S_anchor78','S_anchor81','S_anchorPRO'});

%% Write to Excel
writetable(outTT,'Documents\VARDATA.xlsx','Sheet','S2_Out','WriteMode', 'overwritesheet');
