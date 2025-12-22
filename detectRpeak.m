function [R_locs, R_amp, bp_ecg] = detectRpeak(ecg, Fs, doPlot)
%DETECTRPEAK Detect R-peaks (QRS complexes) in an ECG signal (Pan–Tompkins style).
%
%   [R_locs, R_amp, bp_ecg] = detectRpeak(ecg, Fs)
%   [R_locs, R_amp, bp_ecg] = detectRpeak(ecg, Fs, doPlot)
%
%   Inputs
%     ecg    : ECG signal vector (row/column)
%     Fs     : sampling frequency [Hz]
%     doPlot : true/false, enable diagnostic plots (default false)
%
%   Outputs
%     R_locs : indices of detected R-peaks (samples)
%     R_amp  : amplitudes at R_locs on band-passed ECG
%     bp_ecg : band-passed ECG used for detection

    arguments
        ecg (:,1) double                          % force column, numeric
        Fs (1,1) double {mustBePositive, mustBeFinite}
        doPlot (1,1) logical = false
    end
    
    % --- basic setup ---
    N = numel(ecg);
    t = (0:N-1)'/Fs;
    
    % --- 1) Band-pass filter (QRS enhancement) ---
    bp_ecg = qrs_bandpass(ecg, Fs, doPlot);
    
    if doPlot
        tMax    = min(5, t(end));
        idxPlot = t <= tMax;
    
        figure;
        tiledlayout(2,1);
    
        nexttile;
        plot(t(idxPlot), ecg(idxPlot));
        grid on;
        xlabel('Time (s)'); ylabel('Amplitude');
        title('Original ECG (first 5 s)');
    
        nexttile;
        plot(t(idxPlot), bp_ecg(idxPlot));
        grid on;
        xlabel('Time (s)'); ylabel('Amplitude');
        title('Band-passed ECG (first 5 s)');
    end
    
    % --- 2) Energy envelope: derivative -> square -> moving mean ---
    diff_ecg = derivative(bp_ecg);
    sq_ecg   = diff_ecg.^2;
    
    win_s = 0.15;                                % 150 ms
    win_N = max(1, round(win_s*Fs));
    env_ecg = movmean(sq_ecg, win_N);
    
    if doPlot
        figure;
        tiledlayout(3,1);
        ax = gobjects(3,1);
    
        ax(1) = nexttile;
        plot(t, diff_ecg); grid on;
        xlabel('Time (s)'); ylabel('Amplitude');
        title('Derivative of band-passed ECG');
    
        ax(2) = nexttile;
        plot(t, sq_ecg); grid on;
        xlabel('Time (s)'); ylabel('Amplitude^2');
        title('Squared derivative');
    
        ax(3) = nexttile;
        plot(t, env_ecg); grid on;
        xlabel('Time (s)'); ylabel('Envelope');
        title(sprintf('Moving-window integration (%.0f ms window)', win_s*1e3));
    
        linkaxes(ax,'x');
    end
    
    % --- 3) Peak detection on envelope ---
    init_segment = 1:min(N, round(2*Fs));        % first 2 seconds
    base_thr = 0.3 * max(env_ecg(init_segment) + eps);
    
    minRR_ms   = 200;                             % 200 ms refractory
    minRR_samp = round(minRR_ms*1e-3*Fs);
    
    [pks_env, locs_env] = findpeaks(env_ecg, ...
        'MinPeakHeight', base_thr, ...
        'MinPeakDistance', minRR_samp);
    
    if doPlot
        figure;
        plot(t, env_ecg); hold on;
        plot(t(locs_env), pks_env, 'rv', 'MarkerFaceColor','r');
        grid on;
        xlabel('Time (s)'); ylabel('Envelope');
        title('Envelope with detected peaks (envelope domain)');
        legend('Envelope','Detected peaks');
    end
    
    % --- 4) Refine R locations on band-passed ECG ---
    refine_ms   = 80;
    refine_samp = round(refine_ms*1e-3*Fs);
    
    R_locs = zeros(size(locs_env));
    R_amp  = zeros(size(locs_env));
    
    for k = 1:numel(locs_env)
        ic = locs_env(k);
        i1 = max(1, ic - refine_samp);
        i2 = min(N, ic + refine_samp);
    
        [peak_val, local_idx] = max(bp_ecg(i1:i2));
        R_locs(k) = i1 + local_idx - 1;
        R_amp(k)  = peak_val;
    end
    
    % --- 5) Clean-up ---
    [R_locs, R_amp] = R_cleanup(R_locs, R_amp, Fs);
    
    if doPlot
        figure;
        plot(t, bp_ecg); hold on;
        plot(t(R_locs), bp_ecg(R_locs), 'rv', 'MarkerFaceColor','r');
        grid on;
        xlabel('Time (s)'); ylabel('Amplitude');
        title('Band-passed ECG with detected R-peaks');
        legend('Band-passed ECG','R-peaks');
    end

end

function [R_locs_clean, R_amp_clean] = R_cleanup(R_locs, R_amp, Fs)
%R_CLEANUP  Simple clean-up for R-peaks.
%
%   - Enforce a minimum RR interval
%   - Remove peaks with very small amplitude

    if isempty(R_locs)
        R_locs_clean = R_locs;
        R_amp_clean  = R_amp;
        return;
    end

    % 1) Enforce minimum RR (physiological)
    minRR_ms   = 200;               % 200 ms -> 300 bpm
    minRR_samp = round(minRR_ms*1e-3*Fs);

    keep_idx = true(size(R_locs));
    for k = 2:numel(R_locs)
        if R_locs(k) - R_locs(k-1) < minRR_samp
            % If too close, keep the larger one
            if R_amp(k) > R_amp(k-1)
                keep_idx(k-1) = false;
            else
                keep_idx(k) = false;
            end
        end
    end

    R_locs = R_locs(keep_idx);
    R_amp  = R_amp(keep_idx);

    % 2) Remove very small peaks
    if ~isempty(R_amp)
        med_amp = median(R_amp);
        amp_thr = 0.2 * med_amp;  % keep peaks above 20% of median (tune!)
        valid   = (R_amp >= amp_thr);

        R_locs_clean = R_locs(valid);
        R_amp_clean  = R_amp(valid);
    else
        R_locs_clean = R_locs;
        R_amp_clean  = R_amp;
    end
end