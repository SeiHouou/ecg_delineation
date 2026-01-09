function idx = refineMinBetween(x, Start, End)
%REFINEMINBETWEEN  Find the per-beat minimum index of x within [Start, End].
%
%   idx = refineMinBetween(x, iStart, iEnd)
%
%   Inputs
%   ------
%   x      : signal vector (e.g., raw ECG)
%   Start : start indices (vector, one per beat)
%   End   : end indices   (vector, one per beat)
%
%   Output
%   ------
%   idx    : index of the minimum sample in x(Start:End) for each beat
%            (NaN if the window is invalid or missing)

x = x(:);
Start = Start(:);
End   = End(:);

n = numel(Start);
idx = nan(n,1);

for k = 1:n
    % Skip beats with missing boundaries
    if ~isfinite(Start(k)) || ~isfinite(End(k)), continue; end

    % Clamp boundaries to valid sample range
    a = max(1, round(Start(k)));
    b = min(numel(x), round(End(k)));

    % Invalid/empty window -> leave NaN
    if b <= a, continue; end

    % Find minimum within the beat-specific window
    [~, j] = min(x(a:b));
    idx(k) = a + j - 1;
end
end
