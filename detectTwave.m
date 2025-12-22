function [T_on, T_peak, T_amp, T_off] = detectTwave(ecg, R_locs, S_off, Fs, doPlot, opts)
%DETECTTWAVE  Detect T-wave onset, peak, and offset for each beat.
%
%   [T_on, T_peak, T_amp, T_off] = detectTwave(ecg, R_locs, S_off, Fs)
%   ... = detectTwave(..., doPlot)
%   ... = detectTwave(..., doPlot, opts)
%
%   Inputs
%   ------
%   ecg    : ECG signal (vector)
%   R_locs : R-peak indices (samples)
%   S_off  : QRS offset indices (samples) (same length as R_locs)
%   Fs     : sampling frequency [Hz]
%   doPlot : true/false, diagnostic plots (default false)
%   opts   : struct of parameters (optional)
%
%   Outputs
%   -------
%   T_on   : T-wave onset indices (samples)
%   T_peak : T-wave peak indices (samples)
%   T_amp  : T-wave peak amplitude on ORIGINAL ecg (mV)
%   T_off  : T-wave offset indices (samples)
%
%   Method (simple time-domain heuristic)
%   ------------------------------------
%   1) T-bandpass ECG (default 0.5–10 Hz) to suppress QRS and high-frequency noise.
%   2) Search after S_off in a post-QRS window for the largest |amplitude| -> T_peak.
%   3) From T_peak, walk left/right until both slope (|d/dt|) and amplitude are close
%      to a local baseline estimated from the early ST segment.
%
%   Notes
%   -----
%   This is a lightweight delineation suitable for clean signals. If you have inverted
%   T waves, the |amplitude| criterion still works because we use absolute value.

    arguments
        ecg    {mustBeNumeric, mustBeVector}
        R_locs {mustBeNumeric, mustBeVector}
        S_off  {mustBeNumeric, mustBeVector}
        Fs (1,1) double {mustBePositive, mustBeFinite}
        doPlot (1,1) logical = false

        opts.SearchStart_ms   (1,1) double {mustBeNonnegative} = 60
        opts.SearchEnd_ms     (1,1) double {mustBePositive} = 520
        opts.PreNextR_ms      (1,1) double {mustBeNonnegative} = 120

        opts.OnsetWin_ms      (1,1) double {mustBePositive} = 120
        opts.OffsetWin_ms     (1,1) double {mustBePositive} = 180

        opts.BaseAmpFrac      (1,1) double {mustBePositive} = 0.12
        opts.BaseSlopeFrac    (1,1) double {mustBePositive} = 0.06

        opts.DerivSmooth_ms   (1,1) double {mustBeNonnegative} = 8

        opts.PlotPad_ms       (1,1) double {mustBePositive} = 180
        opts.nShow            (1,1) double {mustBePositive, mustBeInteger} = 3
    end

    % ---- reshape & checks ----
    ecg    = ecg(:);
    R_locs = R_locs(:);
    S_off  = S_off(:);

    N      = numel(ecg);
    nBeats = numel(R_locs);

    if numel(S_off) ~= nBeats
        error('detectTwave:LengthMismatch', ...
            'S_off must have same length as R_locs (got %d vs %d).', numel(S_off), nBeats);
    end

    % ---- outputs ----
    T_on   = NaN(nBeats,1);
    T_peak = NaN(nBeats,1);
    T_amp  = NaN(nBeats,1);
    T_off  = NaN(nBeats,1);

    % ---- preprocess: T-band (reuse p_bandpass with different cutoffs) ----
    bpOpts = struct('LowCut_Hz', 0.5, 'HighCut_Hz', 10, 'Order', 2, 'UseZeroPhase', true);
    t_ecg  = p_bandpass(ecg, Fs, false, bpOpts);

    dt = derivative(t_ecg);
    smoothSamp = max(1, round(opts.DerivSmooth_ms*1e-3*Fs));
    abs_d = movmean(abs(dt), smoothSamp);
    abs_t = abs(t_ecg);

    % ---- ms -> samples ----
    SearchStart_s = round(opts.SearchStart_ms*1e-3*Fs);
    SearchEnd_s   = round(opts.SearchEnd_ms*1e-3*Fs);
    PreNextR_s    = round(opts.PreNextR_ms*1e-3*Fs);

    OnsetWin_s    = round(opts.OnsetWin_ms*1e-3*Fs);
    OffsetWin_s   = round(opts.OffsetWin_ms*1e-3*Fs);

    for k = 1:nBeats
        iR = R_locs(k);

        if iR < 2 || iR > N-1
            continue;
        end

        % Anchor: S_off if available, else use R + 40 ms
        if isfinite(S_off(k)) && S_off(k) >= 1 && S_off(k) <= N
            iSoff = S_off(k);
        else
            iSoff = min(N, iR + round(0.04*Fs));
        end

        i_start = min(N, iSoff + SearchStart_s);
        i_end   = min(N, iR + SearchEnd_s);

        % Prevent spilling into the next beat
        if k < nBeats && isfinite(R_locs(k+1))
            i_end = min(i_end, R_locs(k+1) - PreNextR_s);
        end

        if i_end <= i_start + 10
            continue;
        end

        % --- T peak: max |amplitude| in T-band ---
        [~, rel] = max(abs_t(i_start:i_end));
        iTp = i_start + rel - 1;

        T_peak(k) = iTp;
        T_amp(k)  = ecg(iTp);

        % --- Local baseline: median in early ST segment (S_off..start) ---
        baseSeg = iSoff:min(i_start, N);
        if numel(baseSeg) >= 5
            baseVal = median(t_ecg(baseSeg));
        else
            baseVal = median(t_ecg(max(1, iSoff-20):min(N, iSoff+20)));
        end

        peakAmp = abs(t_ecg(iTp) - baseVal) + eps;
        thrAmp  = opts.BaseAmpFrac * peakAmp;

        localSlopeMax = max(abs_d(i_start:i_end)) + eps;
        thrSlope = opts.BaseSlopeFrac * localSlopeMax;

        % --- T onset: walk left from T_peak ---
        o1 = max(i_start, iTp - OnsetWin_s);
        i = iTp;
        while i > o1
            if abs_d(i) <= thrSlope && abs(t_ecg(i) - baseVal) <= thrAmp
                break;
            end
            i = i - 1;
        end
        T_on(k) = i;

        % --- T offset: walk right from T_peak ---
        o2 = min(i_end, iTp + OffsetWin_s);
        i = iTp;
        while i < o2
            if abs_d(i) <= thrSlope && abs(t_ecg(i) - baseVal) <= thrAmp
                break;
            end
            i = i + 1;
        end
        T_off(k) = i;

        % Sanity: enforce ordering
        if isfinite(T_on(k)) && isfinite(T_off(k)) && T_on(k) >= T_off(k)
            T_on(k)  = NaN;
            T_off(k) = NaN;
        end
    end

    % ---- diagnostic plots ----
    if doPlot
        t = (0:N-1)'/Fs;
        nShow = min(opts.nShow, nBeats);

        for k = 1:nShow
            if ~isfinite(T_peak(k)) || ~isfinite(T_on(k))
                continue;
            end

            padS = round(opts.PlotPad_ms*1e-3*Fs);
            i1 = max(1, S_off(k) - padS);
            i2 = min(N, T_off(k) + padS);

            figure;
            yyaxis left
            plot(t(i1:i2), ecg(i1:i2), 'b', 'DisplayName','ECG'); hold on;
            plot(t(i1:i2), t_ecg(i1:i2), 'k--', 'DisplayName','T-bandpassed ECG');
            grid on;

            plot(t(T_on(k)),   ecg(T_on(k)),   'ko', 'MarkerFaceColor','k', 'DisplayName','T onset');
            plot(t(T_peak(k)), ecg(T_peak(k)), 'ro', 'MarkerFaceColor','r', 'DisplayName','T peak');
            plot(t(T_off(k)),  ecg(T_off(k)),  'go', 'MarkerFaceColor','g', 'DisplayName','T offset');

            xlabel('Time [s]');
            ylabel('ECG [mV]');
            title(sprintf('T-wave delineation – beat %d', k));

            yyaxis right
            plot(t(i1:i2), abs_d(i1:i2), 'm--', 'DisplayName','|dt| (smoothed)');
            ylabel('|dt|');

            legend('show','Location','best');
        end
    end

end
