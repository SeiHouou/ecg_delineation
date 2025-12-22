function out = computeECGIntervals(R_locs, P_on, P_off, Q_on, S_off, T_on, T_off, Fs)
%COMPUTEECGINTERVALS  Compute ECG intervals/segments durations (seconds).
%
%   out = computeECGIntervals(R_locs, P_on, P_off, Q_on, S_off, T_on, T_off, Fs)
%
%   Inputs
%   ------
%   R_locs, P_on, P_off, Q_on, S_off, T_on, T_off : sample indices (vectors)
%   Fs                                            : sampling frequency [Hz]
%
%   Output struct (durations in seconds)
%   -----------------------------------
%   out.RR_s        : RR intervals (nBeats-1 x 1)
%   out.HR_bpm      : Heart rate derived from RR (nBeats-1 x 1)
%   out.PR_s        : PR interval (P_on -> Q_on)
%   out.QRS_s       : QRS duration (Q_on -> S_off)
%   out.QT_s        : QT interval (Q_on -> T_off)
%   out.STint_s     : ST interval (S_off -> T_off)
%   out.PRseg_s     : PR segment (P_off -> Q_on)
%   out.STseg_s     : ST segment (S_off -> T_on)
%
%   Additionally, start/end indices are returned for convenience:
%   out.idx.PR  = [P_on Q_on]
%   out.idx.QRS = [Q_on S_off]
%   out.idx.QT  = [Q_on T_off]
%   out.idx.STi = [S_off T_off]
%   out.idx.PRs = [P_off Q_on]
%   out.idx.STs = [S_off T_on]

    arguments
        R_locs {mustBeNumeric, mustBeVector}
        P_on   {mustBeNumeric, mustBeVector}
        P_off  {mustBeNumeric, mustBeVector}
        Q_on   {mustBeNumeric, mustBeVector}
        S_off  {mustBeNumeric, mustBeVector}
        T_on   {mustBeNumeric, mustBeVector}
        T_off  {mustBeNumeric, mustBeVector}
        Fs (1,1) double {mustBePositive, mustBeFinite}
    end

    % reshape
    R_locs = R_locs(:);
    P_on   = P_on(:);
    P_off  = P_off(:);
    Q_on   = Q_on(:);
    S_off  = S_off(:);
    T_on   = T_on(:);
    T_off  = T_off(:);

    n = numel(R_locs);
    mustHaveSameLength(P_on,n,'P_on');
    mustHaveSameLength(P_off,n,'P_off');
    mustHaveSameLength(Q_on,n,'Q_on');
    mustHaveSameLength(S_off,n,'S_off');
    mustHaveSameLength(T_on,n,'T_on');
    mustHaveSameLength(T_off,n,'T_off');

    % RR and HR
    RR = diff(R_locs) ./ Fs;
    HR = 60 ./ RR;

    % intervals / segments
    PR    = (Q_on  - P_on)  ./ Fs;
    QRS   = (S_off - Q_on)  ./ Fs;
    QT    = (T_off - Q_on)  ./ Fs;
    STint = (T_off - S_off) ./ Fs;
    PRseg = (Q_on  - P_off) ./ Fs;
    STseg = (T_on  - S_off) ./ Fs;

    % sanitize: negatives or non-finite
    PR    = sanitizeDur(PR);
    QRS   = sanitizeDur(QRS);
    QT    = sanitizeDur(QT);
    STint = sanitizeDur(STint);
    PRseg = sanitizeDur(PRseg);
    STseg = sanitizeDur(STseg);

    RR = sanitizeDur(RR);
    HR(~isfinite(RR) | RR<=0) = NaN;

    out = struct();
    out.RR_s    = RR;
    out.HR_bpm  = HR;

    out.PR_s    = PR;
    out.QRS_s   = QRS;
    out.QT_s    = QT;
    out.STint_s = STint;

    out.PRseg_s = PRseg;
    out.STseg_s = STseg;

    out.idx = struct();
    out.idx.PR  = [P_on  Q_on];
    out.idx.QRS = [Q_on  S_off];
    out.idx.QT  = [Q_on  T_off];
    out.idx.STi = [S_off T_off];
    out.idx.PRs = [P_off Q_on];
    out.idx.STs = [S_off T_on];
end

% ----------------- local helpers -----------------
function mustHaveSameLength(x,n,name)
    if numel(x) ~= n
        error('computeECGIntervals:LengthMismatch', '%s must have length %d.', name, n);
    end
end

function d = sanitizeDur(d)
    d(~isfinite(d) | d<=0) = NaN;
end
