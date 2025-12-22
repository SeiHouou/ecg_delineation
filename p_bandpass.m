function y = p_bandpass(x, Fs, doPlot, opts)
%P_BANDPASS  Band-pass filter for P-wave enhancement (~0.5–15 Hz).
%
%   y = p_bandpass(x, Fs)
%   y = p_bandpass(x, Fs, doPlot)
%   y = p_bandpass(x, Fs, doPlot, opts)
%
%   Inputs
%   ------
%   x      : ECG signal (vector)
%   Fs     : sampling frequency [Hz]
%   doPlot : (optional) plot magnitude response (default false)
%   opts   : (optional) parameters:
%            - LowCut_Hz   (default 0.5)
%            - HighCut_Hz  (default 15)
%            - Order       (default 2)  % Butterworth order
%            - UseZeroPhase(default true)  % filtfilt vs filter
%            - PlotMaxHz   (default 60)
%
%   Output
%   ------
%   y : filtered signal (same size as x)
    arguments
        x  {mustBeNumeric, mustBeVector}
        Fs (1,1) double {mustBePositive, mustBeFinite}
        doPlot (1,1) logical = false
        opts struct = struct()
    end

    defaults = struct( ...
        'LowCut_Hz', 0.5, ...
        'HighCut_Hz', 15, ...
        'Order', 2, ...
        'UseZeroPhase', true, ...
        'PlotMaxHz', 60);

    fn = fieldnames(defaults);
    for k = 1:numel(fn)
        if ~isfield(opts, fn{k}) || isempty(opts.(fn{k}))
            opts.(fn{k}) = defaults.(fn{k});
        end
    end
    
    % Keep original shape (row/col)
    x0 = x;
    x  = x(:);
    
    % Validate cutoffs
    if opts.LowCut_Hz >= opts.HighCut_Hz
        error('p_bandpass:BadCutoffs', 'LowCut_Hz must be < HighCut_Hz.');
    end
    if opts.HighCut_Hz >= Fs/2
        error('p_bandpass:BadCutoffs', 'HighCut_Hz must be < Nyquist (Fs/2).');
    end
    
    Wn = [opts.LowCut_Hz opts.HighCut_Hz] / (Fs/2);
    [b,a] = butter(opts.Order, Wn, 'bandpass');
    
    if opts.UseZeroPhase
        y = filtfilt(b,a,x);
    else
        y = filter(b,a,x);
    end
    
    % Restore original shape
    if isrow(x0), y = y.'; end
    
    % Optional plot of filter response
    if doPlot
        Nfft = 4096;
        [H,f] = freqz(b,a,Nfft,Fs);
    
        figure;
        subplot(2,1,1);
        plot(f, 20*log10(abs(H)+eps));
        grid on;
        xlim([0 opts.PlotMaxHz]);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        title(sprintf('P band-pass magnitude (%.2f–%.2f Hz, order %d)', ...
            opts.LowCut_Hz, opts.HighCut_Hz, opts.Order));
    
        subplot(2,1,2);
        plot(f, unwrap(angle(H)));
        grid on;
        xlim([0 opts.PlotMaxHz]);
        xlabel('Frequency (Hz)');
        ylabel('Phase (rad)');
        title('Phase response');
    end

end
