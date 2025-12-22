function y = qrs_bandpass(x, Fs, doPlot)
%QRS_BANDPASS Band-pass filter for QRS enhancement (default 5–25 Hz).
%
%   y = qrs_bandpass(x, Fs)
%   y = qrs_bandpass(x, Fs, doPlot)
%
%   Inputs:
%     x      - ECG signal (vector)
%     Fs     - sampling frequency in Hz
%     doPlot - (optional) true/false to plot frequency response
%
%   Output:
%     y      - filtered ECG (same size as x)

    arguments
        x (:,1) double
        Fs (1,1) double {mustBePositive}
        doPlot (1,1) logical = false
    end
    
    % ---- Parameters (tune if needed) ----
    low_cut  = 5;    % Hz
    high_cut = 25;   % Hz
    order    = 2;    % butter(N, ...) -> bandpass becomes 2N order overall
    
    % ---- Design + apply (zero-phase) ----
    [b, a] = butter(order, [low_cut high_cut]/(Fs/2), "bandpass");
    y = filtfilt(b, a, x);
    
    % ---- Optional plot ----
    if doPlot
        Nfft = 4096;
        [H, f] = freqz(b, a, Nfft, Fs);
    
        figure;
        tiledlayout(2,1);
    
        nexttile;
        plot(f, 20*log10(abs(H)));
        grid on;
        xlim([0 100]);
        xlabel("Frequency (Hz)");
        ylabel("Magnitude (dB)");
        title(sprintf("QRS band-pass magnitude (%.1f–%.1f Hz)", low_cut, high_cut));
    
        nexttile;
        plot(f, unwrap(angle(H)));
        grid on;
        xlim([0 100]);
        xlabel("Frequency (Hz)");
        ylabel("Phase (rad)");
        title("QRS band-pass phase");
    end
end
