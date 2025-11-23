function [r, p, n] = movcorr_edit(x, y, k, varargin)
% MOVCORR  Moving Pearson correlation coefficient
%   r = MOVCORR(x, y, k) returns a column vector r whose t-th entry is
%   the Pearson correlation between the windows {x(t-kB) … x(t+kF)} and
%   {y(t-kB) … y(t+kF)}, where k = [kB kF] or scalar window length.
%   Optional fourth arg is 'includenan' (default) or 'omitnan'.
%
%   [r, p] = MOVCORR(...) also returns p-values; [r, p, n] returns the
%   number of non-NaN observations in each window.
%
%   This version fixes the original File-Exchange bug where n was left
%   empty when missing = 'omitnan' and the window contained no NaNs.

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
[x, y, k, missing, endpoints, tail, nMin, n] = parseInputs(x, y, k, varargin);

isNumericEndpoints = ~ischar(endpoints);
if isNumericEndpoints
    endpointVals = endpoints;
    endpoints    = 'fill';
end

% -------------------------------------------------------------------------
% Rolling sums (using MOVSUM)
% -------------------------------------------------------------------------
movingSum = @(v) movsum(v, k, missing, 'Endpoints', endpoints);

sumX = movingSum(x);
sumY = movingSum(y);

% -------------------------------------------------------------------------
% Pearson r
% -------------------------------------------------------------------------
r = (n .* movingSum(x .* y) - sumX .* sumY) ...
    ./ (sqrt(n .* movingSum(x.^2) - sumX.^2) ...
      .* sqrt(n .* movingSum(y.^2) - sumY.^2));

% Handle values that stray slightly outside ±1 (numerical error)
bad = abs(r) > 1;
r(bad) = r(bad) ./ abs(r(bad));

% Windows with too few data points → NaN
if strcmp(missing, 'omitnan')
    r(n < nMin) = NaN;
end

% -------------------------------------------------------------------------
% p-values
% -------------------------------------------------------------------------
if nargout > 1
    p = computePvalues(r, n, tail);
end

% -------------------------------------------------------------------------
% Restore numeric/logical endpoints, if requested
% -------------------------------------------------------------------------
if isNumericEndpoints
    ids = [1:k(1), numel(r)-k(2)+1:numel(r)];
    r(ids) = endpointVals;
    if nargout > 1, p(ids) = NaN; end
end
end  % main function
% ======================================================================= %
%                          Helper functions                               %
% ======================================================================= %

function [x, y, k, missing, endpoints, tail, nMin, n] = ...
         parseInputs(x, y, k, varArgIn)
% Parse and validate inputs (adapted from original File-Exchange code)

% ---- defaults & validators ---------------------------------------------
VALID_MISS = {'includenan', 'omitnan'};
VALID_ENDP = {'shrink', 'discard', 'fill'};
VALID_TAIL = {'both', 'left', 'right'};

validateCol     = @(v) validateattributes(v, {'numeric'}, {'column'});
validatePosInt  = @(v) validateattributes(v, {'numeric'}, ...
                                          {'scalar','integer','positive'});
validateNonNeg  = @(v) validateattributes(v, {'numeric'}, ...
                                          {'vector','integer','nonnegative'});

% ---- input parser -------------------------------------------------------
ip = inputParser();
ip.StructExpand  = true;
ip.KeepUnmatched = false;

ip.addRequired('x', validateCol);
ip.addRequired('y', validateCol);
ip.addRequired('k', validateNonNeg);
ip.addParameter('Endpoints', VALID_ENDP{1});
ip.addParameter('Tail',      VALID_TAIL{1});
ip.addParameter('MinN',      1, validatePosInt);

% optional fourth positional arg: missing
if isempty(varArgIn) || ~ischar(varArgIn{1}) ...
        || ismember(varArgIn{1}, {'Endpoints','Tail','MinN'})
    missing = VALID_MISS{1};
else
    missing   = varArgIn{1};
    varArgIn(1) = [];
    assert(ismember(missing, VALID_MISS), ...
        'movcorr:MissingInvalid', '''missing'' must be ''omitnan'' or ''includenan''.');
end

ip.parse(x, y, k, varArgIn{:});
endpoints = ip.Results.Endpoints;
tail      = ip.Results.Tail;
nMin      = ip.Results.MinN;

% ---- window length ------------------------------------------------------
if isscalar(k), k = floor(k/2) + [0, -double(mod(k,2)==0)]; end
winSize = sum(k) + 1;

% ---- size checks --------------------------------------------------------
assert(isequal(size(x), size(y)), 'movcorr:XYSize', 'x and y must be same size.');

% ---- initialise n (FIX) -------------------------------------------------
n = winSize;   % ensures n exists even when no NaNs are present

% ---- handle NaNs --------------------------------------------------------
if strcmp(missing, 'omitnan')
    mask = ~(isnan(x) | isnan(y));
    if any(~mask)
        x(~mask) = NaN; y(~mask) = NaN;
        n = movsum(mask, k);
    end
end
end  % parseInputs
% ----------------------------------------------------------------------- %

function p = computePvalues(r, n, tail)
dof = n - 2;
t   = r .* sqrt(dof ./ max(1e-12, 1 - r.^2));  % guard against /0
switch tail
    case 'both',  p = 2 * tcdf(-abs(t), dof);
    case 'left',  p = tcdf(t, dof);
    case 'right', p = tcdf(-t, dof);
end
end
