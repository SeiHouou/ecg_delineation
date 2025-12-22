function [Q_on, Q_peak, R_peak, S_peak, S_off] = delineateQRS(bp_ecg, R_locs, Fs, doPlot, opts)
%DELINEATEQRS  Delineate QRS complex points for each beat.
%
%   [Q_on, Q_peak, R_peak, S_peak, S_off] = delineateQRS(bp_ecg, R_locs, Fs, doPlot)
%   [Q_on, Q_peak, R_peak, S_peak, S_off] = delineateQRS(bp_ecg, R_locs, Fs, doPlot, opts)
%
%   Inputs:
%     bp_ecg  - band-passed ECG (vector)
%     R_locs  - indices of R-peaks (vector)
%     Fs      - sampling frequency [Hz]
%     doPlot  - (optional) true/false, diagnostic plots (default false)
%     opts    - (optional) struct with optional fields:
%               opts.Qon.* overrides Q-on detection params (see below)
%
%   Outputs (all column vectors, length = numel(R_locs)):
%     Q_on    - QRS onset index
%     Q_peak  - Q minimum before R
%     R_peak  - same as R_locs
%     S_peak  - S minimum after R
%     S_off   - QRS offset index

    arguments
        bp_ecg {mustBeNumeric, mustBeVector}
        R_locs {mustBeNumeric, mustBeVector}
        Fs (1,1) double {mustBePositive, mustBeFinite}
        doPlot (1,1) logical = false
        opts struct = struct()
    end

    % ---- reshape / sizes ----
    bp_ecg = bp_ecg(:);
    R_locs = R_locs(:);

    N      = numel(bp_ecg);
    nBeats = numel(R_locs);

    % ---- outputs ----
    Q_on   = NaN(nBeats,1);
    Q_peak = NaN(nBeats,1);
    R_peak = R_locs;
    S_peak = NaN(nBeats,1);
    S_off  = NaN(nBeats,1);

    % store per-beat debug info
    thrSlopeBeat = NaN(nBeats,1);
    iMinBeat     = NaN(nBeats,1);
    QiBeat       = NaN(nBeats,1);

    % ---- parameters (tune here) ----
    qrsHalf_ms    = 80;                       % +/- window around R for Q/S search
    minFracQRStoR = 0.05;                     % reject Q/S if too small vs R

    % ---- Q-on detection parameters (single source of truth) ----
    % Goal: if Q_on is too EARLY (too left), increase onsetSlopeFrac/hold_ms/smooth_ms/baseMadMult
    qon = struct();
    qon.maxPreQ_ms     = 40;    % window length BEFORE Q peak to search for onset
    qon.smooth_ms      = 6;     % smoothing on |dx| used for s
    qon.hold_ms        = 8;     % require s>thr for at least this long
    qon.onsetSlopeFrac = 0.12;  % thr component: fraction of local QRS max slope
    qon.base_ms        = 10;    % baseline slope estimate from start of window
    qon.baseMadMult    = 6;     % baseline threshold = baseMed + baseMadMult*baseMad

    % user overrides: opts.Qon.<field>
    if isfield(opts,'Qon') && ~isempty(opts.Qon)
        f = fieldnames(opts.Qon);
        for ii = 1:numel(f)
            qon.(f{ii}) = opts.Qon.(f{ii});
        end
    end

    % derived sample counts
    qrsHalf_samp = round(qrsHalf_ms*1e-3*Fs);
    maxPreQ_samp = round(qon.maxPreQ_ms*1e-3*Fs);
    smoothN      = max(1, round(qon.smooth_ms*1e-3*Fs));
    holdN        = max(1, round(qon.hold_ms*1e-3*Fs));
    baseN        = max(1, round(qon.base_ms*1e-3*Fs));

    % thresholds for S offset (global, simple)
    dx        = derivative(bp_ecg);

    soff = struct();
    soff.maxPostS_ms = 80;     % don't allow S_off beyond this after S_peak
    soff.slopeFrac   = 0.08;   % fraction of localSlopeMax (bigger => earlier)
    soff.ampFrac     = 0.17;   % fraction of localAmpMax  (bigger => earlier)


    % ---- main loop ----
    for k = 1:nBeats
        iR = R_locs(k);

        % guard: R must be in range
        if iR < 2 || iR > N-1
            continue;
        end

        iL   = max(1, iR - qrsHalf_samp);
        iRgt = min(N, iR + qrsHalf_samp);

        % --- Q peak: min before R ---
        left_seg = bp_ecg(iL:iR);
        [Q_val, idxQ] = min(left_seg);
        Qi = iL + idxQ - 1;

        if abs(Q_val) < minFracQRStoR * abs(bp_ecg(iR))
            Qi = NaN;
        end
        Q_peak(k) = Qi;

        % --- S peak: min after R ---
        right_seg = bp_ecg(iR:iRgt);
        [S_val, idxS] = min(right_seg);
        Si = iR + idxS - 1;

        if abs(S_val) < minFracQRStoR * abs(bp_ecg(iR))
            Si = NaN;
        end
        S_peak(k) = Si;

        % --- Q onset: start of sustained slope rise ---
        if ~isnan(Qi)
            iMin = max(1, Qi - maxPreQ_samp);
            win  = iMin:Qi;  % window where Q_on must lie

            % Local slope scale around the QRS (beat-adaptive)
            qrsWin = max(1, iR - qrsHalf_samp) : min(N, iR + qrsHalf_samp);
            localSlopeMax = max(abs(dx(qrsWin)) + eps);

            % s = smoothed |dx| inside the onset search window
            s = movmean(abs(dx(win)), smoothN);

            % baseline slope estimate from first baseN samples of s
            bIdx    = 1:min(baseN, numel(s));
            baseMed = median(s(bIdx));
            xBase = s(bIdx);
            xBase = xBase(isfinite(xBase));
            if isempty(xBase)
                baseMad = 0;
            else
                m = median(xBase);
                baseMad = median(abs(xBase - m));
            end

            % threshold: be larger of (fraction of QRS max slope) and (baseline+noise margin)
            thr1 = qon.onsetSlopeFrac * localSlopeMax;
            thr2 = baseMed + qon.baseMadMult * baseMad;
            thrSlope = max(thr1, thr2);

            above = (s > thrSlope);

            % sustained-run detection: count trues in trailing window of length holdN
            runCount = movsum(above, [holdN-1 0]);

            jEnd = find(runCount >= holdN, 1, 'first');
            if ~isempty(jEnd)
                jStart = jEnd - holdN + 1;     % start of sustained run
                Qon = win(jStart);
            else
                % fallback: last point still below threshold
                jLastFlat = find(s <= thrSlope, 1, 'last');
                if ~isempty(jLastFlat)
                    Qon = win(jLastFlat);
                else
                    Qon = iMin;
                end
            end

            Q_on(k) = min(Qon, Qi-1);

            % store debug info
            thrSlopeBeat(k) = thrSlope;
            iMinBeat(k)     = iMin;
            QiBeat(k)       = Qi;
        end

        % --- S offset (if S found): beat-adaptive thresholds + max window ---
        if ~isnan(Si)

            maxPostS_samp = round(soff.maxPostS_ms*1e-3*Fs);
            iEnd = min(N, Si + maxPostS_samp);

            % local scales from QRS window around R
            qrsWin = max(1, iR - qrsHalf_samp) : min(N, iR + qrsHalf_samp);
            localSlopeMax = max(abs(dx(qrsWin)) + eps);
            localAmpMax   = max(abs(bp_ecg(qrsWin)) + eps);

            slope_thr_local = soff.slopeFrac * localSlopeMax;
            amp_thr_local   = soff.ampFrac   * localAmpMax;

            i = Si;
            while i < iEnd
                if abs(dx(i)) < slope_thr_local && abs(bp_ecg(i)) < amp_thr_local
                    break;
                end
                i = i + 1;
            end
            S_off(k) = i;
        end

    end

    % ---- plots (optional) ----
    if doPlot
        t = (0:N-1)/Fs;
        nShow = min(3, nBeats);

        for k = 1:nShow
            if isnan(Q_on(k)) || isnan(S_off(k)), continue; end

            pad = round(0.05*Fs);
            i1 = max(1, Q_on(k)  - pad);
            i2 = min(N, S_off(k) + pad);

            figure;

            yyaxis left
            plot(t(i1:i2), bp_ecg(i1:i2), 'b', 'DisplayName','BP ECG');
            hold on; grid on;

            plot(t(Q_on(k)),   bp_ecg(Q_on(k)),   'ko','MarkerFaceColor','k','DisplayName','Q onset');
            plot(t(Q_peak(k)), bp_ecg(Q_peak(k)), 'bo','MarkerFaceColor','b','DisplayName','Q peak');
            plot(t(R_peak(k)), bp_ecg(R_peak(k)), 'ro','MarkerFaceColor','r','DisplayName','R peak');
            plot(t(S_peak(k)), bp_ecg(S_peak(k)), 'mo','MarkerFaceColor','m','DisplayName','S peak');
            plot(t(S_off(k)),  bp_ecg(S_off(k)),  'go','MarkerFaceColor','g','DisplayName','S offset');

            xlabel('Time (s)');
            ylabel('Amplitude');
            title(sprintf('QRS delineation – beat %d', k));

            yyaxis right
            hold on

            % raw |dx| in the shown segment
            dabs = abs(dx(i1:i2));

            % smoothed s over shown segment (same smoothN as algorithm)
            sAll = movmean(dabs, smoothN);

            % mask s outside the Q-on search window (iMin..Qi) so it's clear
            sMask = nan(size(sAll));
            if ~isnan(QiBeat(k)) && ~isnan(iMinBeat(k))
                jA = max(i1, iMinBeat(k));
                jB = min(i2, QiBeat(k));
                if jA <= jB
                    sMask((jA-i1+1):(jB-i1+1)) = sAll((jA-i1+1):(jB-i1+1));
                end
            end

            plot(t(i1:i2), dabs,  'r--', 'DisplayName','|dECG/dt|');
            plot(t(i1:i2), sMask, 'k-',  'LineWidth',1.2, 'DisplayName','s = movmean(|dECG/dt|)');

            thr = thrSlopeBeat(k);
            if isfinite(thr)
                plot([t(i1) t(i2)], [thr thr], 'g-', 'LineWidth',1.5, 'DisplayName','thrSlope');
            end

            ylabel('|dECG/dt|');
            grid on;
            legend('show','Location','best');
        end
    end
end
