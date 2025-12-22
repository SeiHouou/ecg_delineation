function [P_on, P_peak, P_amp, P_off] = detectPwave(ecg, R_locs, Q_on, Fs, doPlot, opts)
%DETECTPWAVE  Detect P-wave onset, peak, and offset for each beat.
%
%   [P_on, P_peak, P_amp]        = detectPwave(ecg, R_locs, Q_on, Fs)
%   [P_on, P_peak, P_amp, P_off] = detectPwave(ecg, R_locs, Q_on, Fs)
%   [...]                        = detectPwave(..., doPlot)
%   [...]                        = detectPwave(..., doPlot, opts)
%
%   Inputs
%   ------
%   ecg    : ECG signal (vector)
%   R_locs : R-peak indices (samples)
%   Q_on   : QRS onset indices (samples) (same length as R_locs)
%   Fs     : sampling frequency [Hz]
%   doPlot : true/false, diagnostic plots (default false)
%   opts   : struct of parameters (optional)
%
%   Outputs
%   -------
%   P_on   : P-wave onset indices (samples)
%   P_peak : P-wave peak indices (samples)
%   P_amp  : P-wave peak amplitudes on ORIGINAL ecg (mV)
%   P_off  : P-wave offset indices (samples)
%
%   Method (simple time-domain heuristic)
%   ------------------------------------
%   1) P-bandpass ECG (default 0.5–15 Hz).
%   2) Use smoothed |d/dt| as a slope indicator.
%   3) P_peak: search a window before Q_on; pick a candidate with small
%      slope and large |amplitude|.
%   4) P_on / P_off: walk left/right from P_peak until both slope and
%      amplitude are close to baseline.
%
%   Notes
%   -----
%   This is a lightweight delineation approach tailored for the project.
%   It should work well for clean, single-lead recordings, but thresholds
%   may require tuning for noisy or atypical morphologies.

    arguments
        ecg    {mustBeNumeric, mustBeVector}
        R_locs {mustBeNumeric, mustBeVector}
        Q_on   {mustBeNumeric, mustBeVector}
        Fs (1,1) double {mustBePositive, mustBeFinite}
        doPlot (1,1) logical = false

        opts.PeakStart_ms     (1,1) double {mustBePositive} = 120
        opts.PeakEnd_ms       (1,1) double {mustBeNonnegative} = 15
        opts.Guard_ms         (1,1) double {mustBeNonnegative} = 10

        opts.OnsetWin_ms      (1,1) double {mustBePositive} = 70
        opts.OffsetWin_ms     (1,1) double {mustBePositive} = 70

        opts.BaseAmpFrac      (1,1) double {mustBePositive} = 0.10
        opts.BaseSlopeFrac    (1,1) double {mustBePositive} = 0.05

        opts.DerivSmooth_ms   (1,1) double {mustBeNonnegative} = 5
        opts.QFallback_ms     (1,1) double {mustBePositive} = 50

        opts.PlotPad_ms       (1,1) double {mustBePositive} = 100
        opts.nShow            (1,1) double {mustBePositive, mustBeInteger} = 3
        opts.BaseWin_ms (1,1) double {mustBePositive} = 20

    end

    % ---- reshape & checks ----
    ecg    = ecg(:);
    R_locs = R_locs(:);
    Q_on   = Q_on(:);

    N      = numel(ecg);
    nBeats = numel(R_locs);

    if numel(Q_on) ~= nBeats
        error('detectPwave:LengthMismatch', ...
            'Q_on must have same length as R_locs (got %d vs %d).', numel(Q_on), nBeats);
    end

    % ---- outputs ----
    P_on   = NaN(nBeats,1);
    P_peak = NaN(nBeats,1);
    P_amp  = NaN(nBeats,1);
    P_off  = NaN(nBeats,1);
    baseValBeat  = NaN(nBeats,1);
    thrAmpBeat   = NaN(nBeats,1);
    thrSlopeBeat = NaN(nBeats,1);

    % ---- preprocess ----
    % Use the dedicated helper in p_bandpass.m (keeps filter parameters centralized).
    p_ecg = p_bandpass(ecg, Fs);

    dp = derivative(p_ecg);
    abs_p = abs(p_ecg);

    smoothSamp = max(1, round(opts.DerivSmooth_ms*1e-3*Fs));
    abs_d = movmean(abs(dp), smoothSamp);

    % ---- ms -> samples ----
    PeakStart_s = round(opts.PeakStart_ms*1e-3*Fs);
    PeakEnd_s   = round(opts.PeakEnd_ms*1e-3*Fs);
    Guard_s     = round(opts.Guard_ms*1e-3*Fs);
    OnsetWin_s  = round(opts.OnsetWin_ms*1e-3*Fs);
    OffsetWin_s = round(opts.OffsetWin_ms*1e-3*Fs);
    BaseWin_s  = round(opts.BaseWin_ms*1e-3*Fs);

    % ---- main loop ----
    for k = 1:nBeats
        iR = R_locs(k);

        if iR < 2 || iR > N-1
            continue;
        end

        % Reference: Q_on if available, otherwise fallback before R.
        if isfinite(Q_on(k)) && Q_on(k) >= 1 && Q_on(k) <= N
            iQ = Q_on(k);
        else
            iQ = max(1, iR - round(opts.QFallback_ms*1e-3*Fs));
        end

        % P-peak search window before QRS
        i_end   = max(1, iQ - max(Guard_s, PeakEnd_s));
        i_start = max(1, iQ - PeakStart_s);

        if i_end <= i_start + 5
            continue;
        end
        
        % ---- local baseline for p_ecg in this beat (use early part of P window) ----
        b1 = i_start;
        b2 = min(i_end, i_start + max(1,BaseWin_s) - 1);
        baseVal = median(p_ecg(b1:b2));

        % deviation from baseline
        dev = abs(p_ecg - baseVal);

        dwin = abs_d(i_start:i_end);

        % Candidate indices: local minima of |dp| (flat selection centered)
        isMin = islocalmin(dwin, 'FlatSelection','center');
        candRel = find(isMin);

        if isempty(candRel)
            [~, j] = min(dwin);
            candRel = j;
        end

        candIdx = i_start + candRel - 1;

        % Choose candidate with largest |amplitude| (avoid baseline traps)
        [~, jBest] = max(dev(candIdx));
        iPp = candIdx(jBest);

        P_peak(k) = iPp;
        P_amp(k)  = ecg(iPp);

        % Local thresholds relative to this beat
        peakDev = dev(iPp) + eps;
        thrAmp  = opts.BaseAmpFrac * peakDev;

        % Slope threshold: relative to the local slope magnitude in search region
        localSlopeMax = max(abs_d(i_start:i_end)) + eps;
        thrSlope = opts.BaseSlopeFrac * localSlopeMax;

        baseValBeat(k)  = baseVal;
        thrAmpBeat(k)   = thrAmp;
        thrSlopeBeat(k) = thrSlope;

        % ---- P onset ----
        o1 = max(i_start, iPp - OnsetWin_s);
        i = iPp;
        while i > o1
            if abs_d(i) <= thrSlope && dev(i) <= thrAmp
                break;
            end
            i = i - 1;
        end
        P_on(k) = i;

        % ---- P offset ----
        f2 = i_end;
        f1 = min(f2, iPp + 1);
        f1 = min(f1, max(1, iPp));

        o2 = min(f2, iPp + OffsetWin_s);
        i = iPp;
        while i < o2
            if abs_d(i) <= thrSlope && dev(i) <= thrAmp
                break;
            end
            i = i + 1;
        end

        % Ensure offset does not run into QRS guard
        P_off(k) = min(i, f2);

        % Sanity: enforce ordering
        if isfinite(P_on(k)) && isfinite(P_off(k)) && P_on(k) >= P_off(k)
            P_on(k)  = NaN;
            P_off(k) = NaN;
        end
    end

    % ---- diagnostic plots ----
    if doPlot
        t = (0:N-1)'/Fs;
        nShow = min(opts.nShow, nBeats);

        for k = 1:nShow
            if ~isfinite(P_peak(k)) || ~isfinite(P_on(k))
                continue;
            end

            padS = round(opts.PlotPad_ms*1e-3*Fs);
            i1 = max(1, P_on(k) - padS);
            i2 = min(N, R_locs(k) + padS);

            figure;
            yyaxis left
            plot(t(i1:i2), ecg(i1:i2), 'b', 'DisplayName','ECG'); hold on;
            plot(t(i1:i2), p_ecg(i1:i2), 'k--', 'DisplayName','P-bandpassed ECG');
            grid on;

            plot(t(P_on(k)),   ecg(P_on(k)),   'ko', 'MarkerFaceColor','k', 'DisplayName','P onset');
            plot(t(P_peak(k)), ecg(P_peak(k)), 'ro', 'MarkerFaceColor','r', 'DisplayName','P peak');
            if isfinite(P_off(k))
                plot(t(P_off(k)), ecg(P_off(k)), 'go', 'MarkerFaceColor','g', 'DisplayName','P offset');
            end

            if isfinite(baseValBeat(k))
                plot([t(i1) t(i2)], [baseValBeat(k) baseValBeat(k)], 'k:', ...
                    'LineWidth', 1.2, 'DisplayName','baseline (p\_ecg)');
            end

            xlabel('Time [s]');
            ylabel('ECG [mV]');
            title(sprintf('P-wave delineation – beat %d', k));

            yyaxis right
            hold on

            % --- recompute dev(t) in the displayed window using stored baseline ---
            baseVal = baseValBeat(k);
            dev_win = abs(p_ecg(i1:i2) - baseVal);

            % Plot slope envelope and its threshold
            plot(t(i1:i2), abs_d(i1:i2)*1e3, 'm--', 'DisplayName','|dp| (smoothed)');
            if isfinite(thrSlopeBeat(k))
                plot([t(i1) t(i2)], [thrSlopeBeat(k) thrSlopeBeat(k)]*1e3, 'm-', ...
                    'LineWidth', 1.5, 'DisplayName','thrSlope');
            end

            % Plot amplitude deviation from baseline and its threshold
            plot(t(i1:i2), dev_win, 'c-', 'LineWidth', 1.2, 'DisplayName','dev = |p\_ecg - baseline|');
            if isfinite(thrAmpBeat(k))
                plot([t(i1) t(i2)], [thrAmpBeat(k) thrAmpBeat(k)], 'c--', ...
                    'LineWidth', 1.5, 'DisplayName','thrAmp');
            end

            ylabel('Slope / deviation');
            grid on


            legend('show','Location','best');
        end
    end

end
