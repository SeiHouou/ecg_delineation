function ax = plotECGWithIntervals(ecg, Fs, titleStr, intervals, markers, varargin)
%PLOTECGWITHINTERVALS  Plot ECG with annotated markers + interval/segment bars.
%
%   ax = plotECGWithIntervals(ecg, Fs, titleStr, intervals, markers)
%   ax = plotECGWithIntervals(..., 'tmin', 0, 'tmax', 10)
%
%   Inputs
%   ------
%   ecg       : ECG vector
%   Fs        : sampling frequency [Hz]
%   titleStr  : title string
%   intervals : output struct from computeECGIntervals(...)
%   markers   : struct with fields (optional): P_on, P_off, P_peak, Q_on, R, S_off, T_on, T_peak, T_off
%
%   Name–Value
%   ----------
%   'tmin', 'tmax'  : time window [s]
%   'Show'          : string array selecting which bars to show
%                    default: ["PR","QRS","QT","STi","PRs","STs"]
%
%   Output
%   ------
%   ax : axes handle

    arguments
        ecg {mustBeNumeric, mustBeVector}
        Fs (1,1) double {mustBePositive, mustBeFinite}
        titleStr {mustBeTextScalar}
        intervals (1,1) struct
        markers (1,1) struct = struct()
    end

    arguments (Repeating)
        varargin
    end

    ecg = ecg(:);
    N   = numel(ecg);
    t   = (0:N-1)'/Fs;

    % --- parse name-value ---
    p = inputParser;
    p.addParameter('tmin', 0, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('tmax', t(end), @(x) isnumeric(x) && isscalar(x));
    p.addParameter('Show', ["PR","QRS","QT","STi","PRs","STs"], @(x) isstring(x) || iscellstr(x));
    p.addParameter('Parent', [], @(h) isempty(h) || isgraphics(h,'axes'));
    
    p.parse(varargin{:});

    tmin = p.Results.tmin;
    tmax = p.Results.tmax;
    show = string(p.Results.Show);
    parentAx = p.Results.Parent;

    % --- plot ECG ---
    if isempty(parentAx)
        figure('Name', char(string(titleStr)));
        ax = axes;
    else
        ax = parentAx;     % use the axes passed from outside
        cla(ax);           % optional: clear existing content
    end

    plot(ax, t, ecg, 'DisplayName','ECG');
    grid(ax,'on'); hold(ax,'on');

    xlim(ax, [tmin tmax]);
    xlabel(ax,'Time [s]');
    ylabel(ax,'ECG Amplitude [mV]');
    title(ax, string(titleStr));

    % --- plot markers (if provided) ---
    plotMarkerField(markers, 'P_on',   'ko', 'P onset');
    plotMarkerField(markers, 'P_peak', 'ro', 'P peak');
    plotMarkerField(markers, 'P_off',  'go', 'P offset');
    plotMarkerField(markers, 'Q_on',   'k^', 'Q onset');
    plotMarkerField(markers, 'R',      'rv', 'R peak');
    plotMarkerField(markers, 'S_off',  'mv', 'S offset');
    plotMarkerField(markers, 'T_on',   'ks', 'T onset');
    plotMarkerField(markers, 'T_peak', 'rs', 'T peak');
    plotMarkerField(markers, 'T_off',  'gs', 'T offset');

    % --- interval bars ---
    barMap = struct(...
        'PR',  struct('name','PR interval', 'idxField','PR'), ...
        'QRS', struct('name','QRS duration','idxField','QRS'), ...
        'QT',  struct('name','QT interval', 'idxField','QT'), ...
        'STi', struct('name','ST interval', 'idxField','STi'), ...
        'PRs', struct('name','PR segment',  'idxField','PRs'), ...
        'STs', struct('name','ST segment',  'idxField','STs'));

    yData = ecg(t>=tmin & t<=tmax);
    if isempty(yData)
        yData = ecg;
    end
    yMin = min(yData);
    yMax = max(yData);
    yRange = max(eps, yMax-yMin);

    baseY = yMin - 0.10*yRange;
    stepY = 0.06*yRange;

    cols = lines(numel(show));
    barHandles = gobjects(numel(show),1);

    for s = 1:numel(show)
        key = show(s);
        if ~isfield(barMap, key)
            continue;
        end

        idx2 = intervals.idx.(barMap.(key).idxField);
        if size(idx2,2) ~= 2
            continue;
        end

        yBar = baseY - (s-1)*stepY;
        first = true;

        for k = 1:size(idx2,1)
            i1 = idx2(k,1);
            i2 = idx2(k,2);
            if ~isfinite(i1) || ~isfinite(i2) || i2 <= i1
                continue;
            end

            tt1 = (i1-1)/Fs;
            tt2 = (i2-1)/Fs;

            % only draw inside viewport
            if tt2 < tmin || tt1 > tmax
                continue;
            end

            tt1 = max(tt1, tmin);
            tt2 = min(tt2, tmax);

            if first
                barHandles(s) = plot(ax, [tt1 tt2], [yBar yBar], '-', ...
                    'LineWidth', 3, 'Color', cols(s,:), 'DisplayName', barMap.(key).name);
                first = false;
            else
                plot(ax, [tt1 tt2], [yBar yBar], '-', 'LineWidth', 3, 'Color', cols(s,:), 'HandleVisibility','off');
            end
        end

        % add a left label once (helps readability)
        if isgraphics(barHandles(s))
            text(ax, tmin, yBar, "  " + barMap.(key).name, 'VerticalAlignment','middle', ...
                'Color', cols(s,:), 'HandleVisibility','off');
        end
    end

    % Expand ylim to include bars
    ylim(ax, [baseY - (numel(show))*stepY - 0.05*yRange, yMax + 0.05*yRange]);

    legend(ax,'show','Location','best');

    % ---- nested helper ----
    function plotMarkerField(S, field, style, name)
        if ~isfield(S, field) || isempty(S.(field)), return; end
        idx = S.(field);
        idx = idx(:);
        idx = idx(isfinite(idx) & idx>=1 & idx<=N);
        if isempty(idx), return; end
        tt = (idx-1)/Fs;
        keep = tt>=tmin & tt<=tmax;
        idx = idx(keep);
        if isempty(idx), return; end
        plot(ax, (idx-1)/Fs, ecg(idx), style, 'MarkerFaceColor','auto', 'DisplayName', name);
    end
end
