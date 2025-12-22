function [resp, env, t] = ecgRespEnvelope(ecg, R_locs, Fs, doPlot, opts)
%ECGRESPENVELOPE  ECG-derived respiration (EDR) + envelope from R-peak amplitudes.
%
%   [resp, env, t] = ecgRespEnvelope(ecg, R_locs, Fs)
%   ... = ecgRespEnvelope(..., doPlot)
%   ... = ecgRespEnvelope(..., doPlot, opts)
%
%   Idea
%   ----
%   Breathing modulates ECG amplitude (chest impedance + heart orientation). A simple
%   EDR (ECG-derived respiration) surrogate is obtained from the beat-to-beat R-peak
%   amplitudes interpolated to a uniform time grid and low-pass filtered.
%
%   Outputs
%   -------
%   resp : respiration surrogate (zero-mean) on a uniform time grid
%   env  : envelope (positive) of resp
%   t    : time vector (seconds), same length as resp/env

    arguments
        ecg    {mustBeNumeric, mustBeVector}
        R_locs {mustBeNumeric, mustBeVector}
        Fs (1,1) double {mustBePositive, mustBeFinite}
        doPlot (1,1) logical = false

        opts.Lowpass_Hz   (1,1) double {mustBePositive} = 0.7
        opts.Order        (1,1) double {mustBePositive, mustBeInteger} = 2
        opts.OutlierMAD   (1,1) double {mustBePositive} = 4
    end

    ecg = ecg(:);
    R_locs = R_locs(:);

    N = numel(ecg);
    t = (0:N-1)'/Fs;

    % --- beat-to-beat amplitude series (EDR proxy) ---
    R_locs = R_locs(isfinite(R_locs) & R_locs>=1 & R_locs<=N);
    R_locs = unique(round(R_locs));

    if numel(R_locs) < 3
        resp = NaN(size(t));
        env  = NaN(size(t));
        return;
    end

    tR   = (R_locs-1)/Fs;
    ampR = ecg(R_locs);

    % --- robust outlier removal on ampR ---
    medA = median(ampR);
    madA = median(abs(ampR - medA)) + eps;
    keep = abs(ampR - medA) <= opts.OutlierMAD * madA;
    tR   = tR(keep);
    ampR = ampR(keep);

    % --- interpolate to uniform grid ---
    ampU = interp1(tR, ampR, t, 'pchip', 'extrap');
    ampU = ampU - mean(ampU, 'omitnan');

    % --- low-pass filter (respiration band) ---
    Wn = opts.Lowpass_Hz / (Fs/2);
    [b,a] = butter(opts.Order, Wn, 'low');
    resp = filtfilt(b,a, ampU);

    % --- envelope (rectify + low-pass again) ---
    envRaw = abs(resp);
    env = filtfilt(b,a, envRaw);

    if doPlot
        figure;
        tiledlayout(2,1);

        nexttile;
        plot(t, ecg); hold on; grid on;
        plot(tR, ampR, 'r.', 'DisplayName','R amplitudes');
        xlabel('Time [s]'); ylabel('ECG [mV]');
        title('ECG with R-peak amplitude samples');
        legend('ECG','R amplitudes');

        nexttile;
        plot(t, resp, 'DisplayName','EDR (resp surrogate)'); hold on; grid on;
        plot(t, env,  'DisplayName','Envelope');
        xlabel('Time [s]'); ylabel('Amplitude (a.u.)');
        title('Respiratory component and envelope (from R-amplitude modulation)');
        legend('show','Location','best');
    end
end
