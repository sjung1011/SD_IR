% Clean and Create Data for inflation-relateds regimes: 
% Good and Bad inflations (S1)
% FRED data
    % Monthly: PCE(chain-type, 2017=100, (1959 Jan-2008 Jan)), CPI (1957 Jan-2008 Jan) 
    % Quarterly: real nondurables and real services quantity indices

clc
clear all
close all

    
%% Load the spreadsheet
Tm = readtable('VARDATA.xlsx','Sheet','S1_M');      
Tq = readtable('VARDATA.xlsx','Sheet','S1_Q');
%% back-splice PCE with CPI during the missing period 
pce  = Tm.PCEPI;                % NaN before 1959-01
cpi  = Tm.CPIAUCSL;             % CPI index, 1982-84 = 100, always available
dln_cpi = [NaN; diff(log(cpi))];       % Δln CPI_t   (first obs = NaN)

% --- Initialise synthetic log-PCEPI ---------------------------------------
ln_pce_syn        = NaN(size(pce));           % start with actual PCEPI (NaN early)
idx = find(Tm.observation_date == '01-Jan-1959');
ln_pce_syn(idx) = log(pce(idx));       % ***anchor at the actual value***

% --- Back-cast ------------------------------------------------------------
for k = idx:-1:2
    if isnan(ln_pce_syn(k-1))                % only fill missing slots
        ln_pce_syn(k-1) = ln_pce_syn(k) - dln_cpi(k);
    end
end

% --- Combine real & synthetic --------------------------------------------
pce_backcast         = exp(ln_pce_syn);     % level series 1957-today
pce_complete         = pce;                         % copy original series
need_fill            = isnan(pce_complete);         % true before 1959-01
pce_complete(need_fill) = pce_backcast(need_fill);    % splice in back-cast
Tm.PCEPI_filled = pce_complete;

%% Convert to a timetable
%   Assume the sheet contains a column called 'DATE' with YYYY-MM-DD strings
Tm_tt = table2timetable(Tm,'RowTimes','observation_date'); % monthly
Tq_tt = table2timetable(Tq,'RowTimes','observation_date'); % quaterly
Tq_monthly = retime(Tq_tt,'monthly','previous');    % step-interpolate; ow 'linear' for a smooth interpolation.
% Combine all to one timetable 
bigTT = synchronize(Tm_tt, Tq_monthly, 'union','previous'); % fill missing values with previous observation

%% Construct demand and inflation growth rates
% Demand 
    % Log-levels, index
bigTT.lnC_id = log(bigTT.RNDURID) + log(bigTT.RSERVID);      
    % Per capita
bigTT.lnC_pc = log(bigTT.RNDURPC) + log(bigTT.RSERVPC);
    % Year-over-year (YoY) growth rates  
bigTT.dC_id  = 100*(bigTT.lnC_id - lagmatrix(bigTT.lnC_id,12)); % or MoM
bigTT.dC_pc  = 100*(bigTT.lnC_pc - lagmatrix(bigTT.lnC_pc,12));

% Inflation (PCEPI YoY)
bigTT.dP  = 100*(log(bigTT.PCEPI_filled) - lagmatrix(log(bigTT.PCEPI_filled),12));

%% Compute the rolling correlation
% Extract the sub-sample you care about
t0 = datetime(1962,7,1);        % start of regime indicator
t1 = datetime(2008,6,1);        % end
subTT = bigTT(t0:t1, :);
% Rolling 60-month Pearson correlation
% movcorr(x,y,[w-1 0]) uses data from t-w+1 through t, i.e. a trailing window.
% index
w  = 60;                        % window length 
subTT.rCorrid_60 = movcorr_edit(subTT.dP, subTT.dC_id, [w-1 0], 'omitnan');
w = 30;
subTT.rCorrid_30 = movcorr_edit(subTT.dP, subTT.dC_id, [w-1 0], 'omitnan');
% per capita
w  = 60;                        % window length 
subTT.rCorrpc_60 = movcorr_edit(subTT.dP, subTT.dC_pc, [w-1 0], 'omitnan');
w = 30;
subTT.rCorrpc_30 = movcorr_edit(subTT.dP, subTT.dC_pc, [w-1 0], 'omitnan');

% Create and define the binary state S
    % Good (demand-driven) inflation (1): positive co-movement of demand and prices.
    % Bad (supply-driven) inflation (0): negative co-movement.
    % zero threshold
subTT.S1_60_id = double(subTT.rCorrid_60 > 0);  
subTT.S1_30_id = double(subTT.rCorrid_30 > 0);
subTT.S1_60_pc = double(subTT.rCorrpc_60 > 0);  
subTT.S1_30_pc = double(subTT.rCorrpc_30 > 0);

%% Count each regime
% --- 60-month indicator ---------------------------------------------------
n1_60 = nnz(subTT.S1_60_id == 1);                       % regime 1
n0_60 = nnz(subTT.S1_60_id == 0);                       % regime 0

% --- 30-month indicator ---------------------------------------------------
n1_30 = nnz(subTT.S1_30_id == 1);
n0_30 = nnz(subTT.S1_30_id == 0);

fprintf('60-month window:  Regime 1 = %d   Regime 0 = %d\n', n1_60, n0_60);
fprintf('30-month window:  Regime 1 = %d   Regime 0 = %d\n', n1_30, n0_30);
% 60-month window:  Regime 1 = 185   Regime 0 = 367
% 30-month window:  Regime 1 = 237   Regime 0 = 315

% --- 60-month indicator ---------------------------------------------------
n1_60 = nnz(subTT.S1_60_pc == 1);                       % regime 1
n0_60 = nnz(subTT.S1_60_pc == 0);                       % regime 0

% --- 30-month indicator ---------------------------------------------------
n1_30 = nnz(subTT.S1_30_pc == 1);
n0_30 = nnz(subTT.S1_30_pc == 0);

fprintf('60-month window:  Regime 1 = %d   Regime 0 = %d\n', n1_60, n0_60);
fprintf('30-month window:  Regime 1 = %d   Regime 0 = %d\n', n1_30, n0_30);
% 60-month window:  Regime 1 = 191   Regime 0 = 361
% 30-month window:  Regime 1 = 239   Regime 0 = 313

%% Save and export 
outTT = subTT(:,{'PCEPI_filled','dP','dC_id','dC_pc','rCorrid_30','rCorrid_60', ...
    'rCorrpc_30','rCorrpc_60',...
    'S1_60_id','S1_60_pc','S1_30_id','S1_30_pc'});

%% Write to Excel
writetimetable(outTT,'Documents\VARDATA.xlsx','Sheet','S1_Out','WriteMode', 'overwritesheet');
