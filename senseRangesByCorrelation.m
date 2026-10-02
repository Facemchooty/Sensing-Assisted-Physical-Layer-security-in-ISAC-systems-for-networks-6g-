function rngEst = senseRangesByCorrelation(txTime, bsPos, obj, fc, Fs)
%SENSERANGESBYCORRELATION Rough monostatic range estimation using correlation
%
% txTime : transmitted time-domain samples (column)
% bsPos  : 2x1 base station position [x;y]
% obj    : struct array with fields .pos (2x1) and optionally .rcs
% fc     : carrier frequency (Hz) [not used directly here, kept for consistency]
% Fs     : sampling rate (Hz)
%
% Output:
% rngEst : 1xNo estimated ranges (meters)

    c = 3e8;
    x = txTime(:);
    No = numel(obj);
    rngEst = nan(1, No);

    for i = 1:No
        d = norm(obj(i).pos - bsPos);
        tau = 2*d/c;                    % two-way delay
        n0 = round(tau * Fs);

        if n0 >= numel(x)
            continue;
        end

        rcs = 1;
        if isfield(obj(i),'rcs')
            rcs = max(obj(i).rcs, 0);
        end

        echo = [zeros(n0,1); x(1:end-n0)] * (0.1 * sqrt(rcs));

        r = xcorr(echo, x);
        [~,ix] = max(abs(r));
        lag = ix - numel(x);

        tauHat = abs(lag) / Fs;
        rngEst(i) = (tauHat * c) / 2;
    end
end