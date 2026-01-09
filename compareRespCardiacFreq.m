function out = compareRespCardiacFreq(resp, Fs, R_locs, doPlot, opts)
%COMPARERESPCARDIACFREQ  Estimate respiratory and cardiac frequencies & ranges.
%
%   out = compareRespCardiacFreq(resp, Fs, R_locs)
%   out = compareRespCardiacFreq(..., doPlot)
%   out = compareRespCardiacFreq(..., doPlot, opts)
%
%   Inputs
%   ------
%   resp   : respiration surrogate signal (uniformly sampled)
%   Fs     : sampling frequency [Hz]
%   R_locs : R-peak indices for cardiac rate estimation
%
%   Output fields
%   -------------
%   out.resp_Hz_dominant  : dominant respiration frequency (Hz) via Welch PSD
%   out.resp_Hz_range     : [min max] breathing frequency from peak-to-peak times
%   out.card_Hz_mean      : mean cardiac frequency (Hz) from RR intervals
%   out.card_Hz_range     : [min max] cardiac frequency (Hz) from RR intervals
%
%   Notes
%   -----
%   Typical ranges:
%     respiration: ~0.1–0.5 Hz (6–30 breaths/min)
%     cardiac:     ~0.8–2 Hz   (48–120 bpm)

    arguments
        resp  {mustBeNumeric, mustBeVector}
        Fs (1,1) double {mustBePositive, mustBeFinite}
        R_locs {mustBeNumeric, mustBeVector}
        doPlot (1,1) logical = false

        opts.RespBand_Hz (1,2) double {mustBePositive} = [0.05 0.8]
        opts.MinBreath_s (1,1) double {mustBePositive} = 1.2
        opts.MaxBreath_s (1,1) double {mustBePositive} = 12
        opts.RespFs_Hz   (1,1) double {mustBePositive} = 25
    end

    resp = resp(:);
    R_locs = R_locs(:);

    % --- Cardiac frequency from RR intervals --- (use original Fs)
    RR_s = diff(R_locs) ./ Fs;
    RR_s = RR_s(isfinite(RR_s) & RR_s > 0);

    card_Hz = 1 ./ RR_s;
    out.card_Hz_mean  = mean(card_Hz, 'omitnan');
    out.card_Hz_range = [min(card_Hz,[],'omitnan') max(card_Hz,[],'omitnan')];

    % --- Respiration (optional sub-sampling for low-frequency analysis) ---
    Fs_resp = Fs;
    if opts.RespFs_Hz < Fs
        ds = max(1, round(Fs / opts.RespFs_Hz));
        resp = resp(1:ds:end);
        Fs_resp = Fs / ds;
    end

    % --- Respiration dominant frequency from PSD ---
    n = numel(resp);
    win = max(256, round(8*Fs_resp));
    win = min(win, n);
    if mod(win,2)==1, win = win-1; end

    [Pxx,f] = pwelch(resp, hamming(win), round(0.5*win), [], Fs_resp);
    band = f>=opts.RespBand_Hz(1) & f<=opts.RespBand_Hz(2);
    if any(band)
        [~,iMax] = max(Pxx(band));
        fBand = f(band);
        out.resp_Hz_dominant = fBand(iMax);
    else
        out.resp_Hz_dominant = NaN;
    end

    % --- Respiration range from peak-to-peak times (time domain) ---
    minDist = round(opts.MinBreath_s * Fs_resp);
    [~,locs] = findpeaks(resp, 'MinPeakDistance', minDist);

    if numel(locs) >= 3
        T = diff(locs) ./ Fs_resp;
        T = T(T>=opts.MinBreath_s & T<=opts.MaxBreath_s);
        resp_Hz = 1 ./ T;
        out.resp_Hz_range = [min(resp_Hz,[],'omitnan') max(resp_Hz,[],'omitnan')];
    else
        out.resp_Hz_range = [NaN NaN];
    end

    if doPlot
        figure;
        plot(f, Pxx); grid on;
        xlim([0 5]);
        xlabel('Frequency [Hz]'); ylabel('PSD');
        title('Respiration surrogate PSD (Welch)');
        hold on;
        if isfinite(out.resp_Hz_dominant)
            xline(out.resp_Hz_dominant, '--', 'Dominant f_{resp}', 'LabelVerticalAlignment','bottom');
        end
    end
end
