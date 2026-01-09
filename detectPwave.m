function [P_on, P_peak, P_amp, P_off] = detectPwave(ecg, R_locs, Q_on, Fs, doPlot, opts)
%DETECTPWAVE  Detect P-wave onset, peak, and offset for each beat.
%
% Algorithm (per beat)
% 1) Filter ECG to emphasize P-waves: p_ecg = p_bandpass(ecg)
% 2) Build slope envelope: abs_d = movmean(abs(derivative(p_ecg)))
% 3) Define a P-search window before QRS (anchored to Q_on, else fallback before R)
% 4) Estimate local baseline from early part of the P-search window
%    dev(t) = |p_ecg(t) - baseline|
% 5) P_peak: choose candidate points with small slope (local minima of abs_d)
%    and pick the one with the largest dev (strongest atrial deflection)
% 6) P_on: walk LEFT from P_peak until (abs_d <= thrSlope_on) AND (dev <= thrAmp_on)
% 7) P_off: walk RIGHT from P_peak until (abs_d <= thrSlope_off) AND (dev <= thrAmp_off)

    arguments
        ecg    {mustBeNumeric, mustBeVector}
        R_locs {mustBeNumeric, mustBeVector}
        Q_on   {mustBeNumeric, mustBeVector}
        Fs (1,1) double {mustBePositive, mustBeFinite}
        doPlot (1,1) logical = false

        % --- windows (ms) ---
        opts.PeakStart_ms (1,1) double {mustBePositive}    = 120  % search start before Q
        opts.PeakEnd_ms   (1,1) double {mustBeNonnegative} = 15   % search end before Q
        opts.Guard_ms     (1,1) double {mustBeNonnegative} = 10   % keep away from QRS

        opts.OnsetWin_ms  (1,1) double {mustBePositive}    = 70   % max onset lookback from P_peak
        opts.OffsetWin_ms (1,1) double {mustBePositive}    = 70   % max offset lookahead from P_peak

        % --- shared threshold fractions ---
        opts.BaseAmpFrac   (1,1) double {mustBePositive} = 0.10
        opts.BaseSlopeFrac (1,1) double {mustBePositive} = 0.05

        % --- optional separate tuning (NaN => use shared) ---
        opts.BaseAmpFrac_on    (1,1) double = NaN
        opts.BaseSlopeFrac_on  (1,1) double = NaN
        opts.BaseAmpFrac_off   (1,1) double = NaN
        opts.BaseSlopeFrac_off (1,1) double = NaN

        % --- preprocessing ---
        opts.DerivSmooth_ms (1,1) double {mustBeNonnegative} = 5
        opts.BaseWin_ms     (1,1) double {mustBePositive}    = 20
        opts.QFallback_ms   (1,1) double {mustBePositive}    = 50

        % --- plotting ---
        opts.PlotPad_ms     (1,1) double {mustBePositive} = 100
        opts.nShow          (1,1) double {mustBePositive, mustBeInteger} = 3
        opts.PlotSlopeScale (1,1) double {mustBePositive} = 1e3 % purely for visibility
    end

    % -------------------------
    % Input shaping & checks
    % -------------------------
    ecg    = ecg(:);
    R_locs = R_locs(:);
    Q_on   = Q_on(:);

    N      = numel(ecg);
    nBeats = numel(R_locs);

    if numel(Q_on) ~= nBeats
        error('detectPwave:LengthMismatch', ...
            'Q_on must have same length as R_locs (got %d vs %d).', numel(Q_on), nBeats);
    end

    % -------------------------
    % Preprocess
    % -------------------------
    p_ecg = p_bandpass(ecg, Fs);
    dp    = derivative(p_ecg);

    smoothN = max(1, round(opts.DerivSmooth_ms*1e-3*Fs));
    abs_d   = movmean(abs(dp), smoothN);   % slope envelope

    % -------------------------
    % ms -> samples
    % -------------------------
    PeakStart_s = round(opts.PeakStart_ms*1e-3*Fs);
    PeakEnd_s   = round(opts.PeakEnd_ms*1e-3*Fs);
    Guard_s     = round(opts.Guard_ms*1e-3*Fs);
    OnsetWin_s  = round(opts.OnsetWin_ms*1e-3*Fs);
    OffsetWin_s = round(opts.OffsetWin_ms*1e-3*Fs);
    BaseWin_s   = round(opts.BaseWin_ms*1e-3*Fs);
    Qfb_s       = round(opts.QFallback_ms*1e-3*Fs);

    % -------------------------
    % Outputs
    % -------------------------
    P_on   = NaN(nBeats,1);
    P_peak = NaN(nBeats,1);
    P_amp  = NaN(nBeats,1);
    P_off  = NaN(nBeats,1);

    % Debug storage (for plots)
    dbg.baseVal   = NaN(nBeats,1);
    dbg.i_start   = NaN(nBeats,1);
    dbg.i_end     = NaN(nBeats,1);
    dbg.o1        = NaN(nBeats,1);
    dbg.o2        = NaN(nBeats,1);
    dbg.thrAmp_on = NaN(nBeats,1);
    dbg.thrSlp_on = NaN(nBeats,1);
    dbg.thrAmp_off= NaN(nBeats,1);
    dbg.thrSlp_off= NaN(nBeats,1);

    % -------------------------
    % Main loop (per beat)
    % -------------------------
    for k = 1:nBeats
        iR = R_locs(k);
        if iR < 2 || iR > N-1
            continue;
        end

        % Anchor window to Q_on (preferred), else fallback before R
        if isfinite(Q_on(k)) && Q_on(k) >= 1 && Q_on(k) <= N
            iQ = Q_on(k);
        else
            iQ = max(1, iR - Qfb_s);
        end

        % P-peak search window BEFORE QRS
        i_end   = max(1, iQ - max(Guard_s, PeakEnd_s));
        i_start = max(1, iQ - PeakStart_s);

        if i_end <= i_start + 5
            continue;
        end

        % Baseline from early part of P window (PR segment)
        b1 = i_start;
        b2 = min(i_end, i_start + max(1,BaseWin_s) - 1);
        baseVal = median(p_ecg(b1:b2));

        % Deviation from baseline (used for "back to baseline" logic)
        dev = abs(p_ecg - baseVal);

        % --- P_peak candidates: local minima of slope envelope in the search window ---
        dwin = abs_d(i_start:i_end);
        isMin = islocalmin(dwin, 'FlatSelection','center');
        candRel = find(isMin);

        if isempty(candRel)
            [~, j] = min(dwin);
            candRel = j;
        end
        candIdx = i_start + candRel - 1;

        % Choose the candidate with strongest deviation from baseline
        [~, jBest] = max(dev(candIdx));
        iPp = candIdx(jBest);

        P_peak(k) = iPp;
        P_amp(k)  = ecg(iPp);

        % --- Thresholds (separated onset vs offset) ---
        peakDev = dev(iPp) + eps;
        localSlopeMax = max(abs_d(i_start:i_end)) + eps;

        [ampFrac_on,  slpFrac_on ] = pickFracs(opts.BaseAmpFrac, opts.BaseSlopeFrac, ...
                                               opts.BaseAmpFrac_on, opts.BaseSlopeFrac_on);

        [ampFrac_off, slpFrac_off] = pickFracs(opts.BaseAmpFrac, opts.BaseSlopeFrac, ...
                                               opts.BaseAmpFrac_off, opts.BaseSlopeFrac_off);

        thrAmp_on   = ampFrac_on  * peakDev;
        thrSlp_on   = slpFrac_on  * localSlopeMax;
        thrAmp_off  = ampFrac_off * peakDev;
        thrSlp_off  = slpFrac_off * localSlopeMax;

        % --- P_on (walk left) ---
        o1 = max(i_start, iPp - OnsetWin_s);
        i  = iPp;
        while i > o1
            if abs_d(i) <= thrSlp_on && dev(i) <= thrAmp_on
                break;
            end
            i = i - 1;
        end
        P_on(k) = i;

        % --- P_off (walk right) ---
        o2 = min(i_end, iPp + OffsetWin_s);
        % i  = iPp;
        % while i < o2
        %     if abs_d(i) <= thrSlp_off && dev(i) <= thrAmp_off
        %         break;
        %     end
        %     i = i + 1;
        %  end
        % P_off(k) = min(i, i_end);
        idx = iPp:o2;
        validOff = (abs_d(idx) <= thrSlp_off) & (dev(idx) <= thrAmp_off);
        if any(validOff)
            P_off(k) = idx(find(validOff, 1, 'last'));
        else
            P_off(k) = o2;
        end
        % Sanity
        if isfinite(P_on(k)) && isfinite(P_off(k)) && P_on(k) >= seeFinite(P_off(k))
            P_on(k)  = NaN;
            P_off(k) = NaN;
        end

        % store debug
        dbg.baseVal(k)   = baseVal;
        dbg.i_start(k)   = i_start;
        dbg.i_end(k)     = i_end;
        dbg.o1(k)        = o1;
        dbg.o2(k)        = o2;
        dbg.thrAmp_on(k) = thrAmp_on;
        dbg.thrSlp_on(k) = thrSlp_on;
        dbg.thrAmp_off(k)= thrAmp_off;
        dbg.thrSlp_off(k)= thrSlp_off;
    end

    % -------------------------
    % Debug plots
    % -------------------------
    if doPlot
        t = (0:N-1)'/Fs;
        nShow = min(opts.nShow, nBeats);
        padS  = round(opts.PlotPad_ms*1e-3*Fs);
        scl   = opts.PlotSlopeScale;

        for k = 1:nShow
            if ~isfinite(P_peak(k)) || ~isfinite(P_on(k)) || ~isfinite(P_off(k))
                continue;
            end

            % --- ONSET figure (zoom around onset->peak) ---
            i1 = max(1, dbg.o1(k) - padS);
            i2 = min(N, P_peak(k) + padS);
            plotPdebug(t, ecg, p_ecg, abs_d, dbg, P_on, P_peak, P_off, k, i1, i2, scl, "onset");

            % --- OFFSET figure (zoom around peak->offset) ---
            i1 = max(1, P_peak(k) - padS);
            i2 = min(N, dbg.i_end(k) + padS);
            plotPdebug(t, ecg, p_ecg, abs_d, dbg, P_on, P_peak, P_off, k, i1, i2, scl, "offset");
        end
    end
end

% =====================================================================
% Helpers (local functions)
% =====================================================================

function [ampFrac, slopeFrac] = pickFracs(ampDefault, slopeDefault, ampOverride, slopeOverride)
% Use overrides if provided (not NaN), else use defaults.
    ampFrac   = ampDefault;
    slopeFrac = slopeDefault;
    if ~isnan(ampOverride),   ampFrac   = ampOverride;   end
    if ~isnan(slopeOverride), slopeFrac = slopeOverride; end
    if ampFrac <= 0 || slopeFrac <= 0
        error('detectPwave:BadFrac', 'Amplitude/slope fractions must be > 0.');
    end
end

function plotPdebug(t, ecg, p_ecg, abs_d, dbg, P_on, P_peak, P_off, k, i1, i2, scl, modeStr)
% Two separate debug figures:
%  - modeStr="onset": show *_on thresholds
%  - modeStr="offset": show *_off thresholds

    baseVal = dbg.baseVal(k);
    dev_win = abs(p_ecg(i1:i2) - baseVal);

    if modeStr == "onset"
        thrAmp = dbg.thrAmp_on(k);
        thrSlp = dbg.thrSlp_on(k);
        titleStr = sprintf('P-wave ONSET debug – beat %d', k);
        xMark = P_on(k);
        xLabel = 'P onset';
        xLimLine = dbg.o1(k);
        xLimName = 'o1 (onset limit)';
    else
        thrAmp = dbg.thrAmp_off(k);
        thrSlp = dbg.thrSlp_off(k);
        titleStr = sprintf('P-wave OFFSET debug – beat %d', k);
        xMark = P_off(k);
        xLabel = 'P offset';
        xLimLine = dbg.o2(k);
        xLimName = 'o2 (offset limit)';
    end

    figure;
    yyaxis left
    plot(t(i1:i2), ecg(i1:i2), 'b', 'DisplayName','ECG'); hold on; grid on;
    plot(t(i1:i2), p_ecg(i1:i2), 'k--', 'DisplayName','p\_ecg (P-bandpass)');
    plot([t(i1) t(i2)], [baseVal baseVal], 'k:', 'LineWidth',1.2, 'DisplayName','baseline (p\_ecg)');

    plot(t(P_peak(k)), ecg(P_peak(k)), 'ro', 'MarkerFaceColor','r', 'DisplayName','P peak');
    plot(t(P_on(k)),   ecg(P_on(k)),   'ko', 'MarkerFaceColor','k', 'DisplayName','P onset');
    plot(t(P_off(k)),  ecg(P_off(k)),  'go', 'MarkerFaceColor','g', 'DisplayName','P offset');

    if isfinite(xLimLine)
        xline(t(xLimLine), 'k--', 'DisplayName', xLimName);
    end

    xlabel('Time [s]');
    ylabel('ECG [mV]');
    title(titleStr);

    yyaxis right
    plot(t(i1:i2), abs_d(i1:i2)*scl, 'm--', 'DisplayName',sprintf('|dp| (smoothed) x%.0g', scl));
    plot([t(i1) t(i2)], [thrSlp thrSlp]*scl, 'm-', 'LineWidth',1.5, 'DisplayName','thrSlope');
    plot(t(i1:i2), dev_win, 'c-', 'LineWidth',1.2, 'DisplayName','dev = |p\_ecg - baseline|');
    plot([t(i1) t(i2)], [thrAmp thrAmp], 'c--', 'LineWidth',1.5, 'DisplayName','thrAmp');

    ylabel('Slope / deviation');
    legend('show','Location','best');
end

function y = seeFinite(x)
% tiny helper: keeps the ordering sanity-check robust if NaN
    if isfinite(x), y = x; else, y = inf; end
end
