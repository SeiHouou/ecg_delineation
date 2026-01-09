function [resp, env, t] = ecgRespLowpass(ecg, Fs, doPlot, opts)
%ECGRESPLOWPASS  Estimate respiratory component via low-pass ECG filtering.
%
%   [resp, env, t] = ecgRespLowpass(ecg, Fs)
%   ... = ecgRespLowpass(ecg, Fs, doPlot)
%   ... = ecgRespLowpass(ecg, Fs, doPlot, opts)
%
%   Inputs
%   ------
%   ecg    : ECG signal (vector)
%   Fs     : sampling frequency [Hz]
%   doPlot : true/false, diagnostic plots (default false)
%   opts   : struct of parameters (optional)
%
%   Outputs
%   -------
%   resp : respiration surrogate (low-pass ECG)
%   env  : envelope of resp (positive, smoothed)
%   t    : time vector (seconds), same length as resp/env
%
%   Notes
%   -----
%   This is a simple alternative to R-peak amplitude modulation. It assumes
%   the respiratory component lies below the cutoff frequency.

    arguments
        ecg    {mustBeNumeric, mustBeVector}
        Fs (1,1) double {mustBePositive, mustBeFinite}
        doPlot (1,1) logical = false
    
        opts.Lowpass_Hz (1,1) double {mustBePositive} = 0.7
        opts.Order      (1,1) double {mustBePositive, mustBeInteger} = 2
    end
    
    ecg = ecg(:);
    N = numel(ecg);
    t = (0:N-1)'/Fs;
    
    % --- low-pass filter (respiration band) ---
    Wn = opts.Lowpass_Hz / (Fs/2);
    [b,a] = butter(opts.Order, Wn, 'low');
    resp = filtfilt(b,a, ecg);
    resp = resp - mean(resp, 'omitnan');
    
    % --- envelope (rectify + low-pass again) ---
    envRaw = abs(resp);
    env = filtfilt(b,a, envRaw);
    
    if doPlot
        figure;
        tiledlayout(2,1);
    
        nexttile;
        plot(t, ecg, 'DisplayName','ECG'); grid on;
        xlabel('Time [s]'); ylabel('ECG [mV]');
        title('ECG (original)');
    
        nexttile;
        plot(t, resp, 'DisplayName','Resp (low-pass)'); hold on; grid on;
        plot(t, env,  'DisplayName','Envelope');
        xlabel('Time [s]'); ylabel('Amplitude (a.u.)');
        title('Respiratory component and envelope (low-pass)');
        legend('show','Location','best');
    end
    
end
