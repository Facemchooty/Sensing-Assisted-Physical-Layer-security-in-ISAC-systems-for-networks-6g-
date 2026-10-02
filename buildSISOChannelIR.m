function h = buildSISOChannelIR(txPos, rxPos, obj, fc, Fs)
%BUILDSISOCHANNELIR Simple SISO multipath impulse response (2D)
% LoS + a few reflected paths via moving objects (tx->obj->rx)
%
% txPos, rxPos: 2x1 vectors [x;y]
% obj: struct array with fields .pos (2x1), .rcs
% fc: carrier frequency (Hz)
% Fs: sampling rate (Hz)
%
% Output:
% h: column FIR taps (complex), delays quantized to samples

    c = 3e8;

    % ---- LoS ----
    d0   = norm(rxPos - txPos);
    tau0 = d0 / c;
    a0   = friisAmp(d0, fc) * exp(1j*2*pi*rand);
    n0   = round(tau0 * Fs);

    taps  = n0;
    gains = a0;

    % ---- Reflections via objects ----
    K = min(3, numel(obj));  % keep small for speed
    for i = 1:K
        d1  = norm(obj(i).pos - txPos);
        d2  = norm(rxPos - obj(i).pos);
        d   = d1 + d2;
        tau = d / c;
        n   = round(tau * Fs);

        % smaller gain for reflection
        a = 0.25 * friisAmp(d, fc) * exp(1j*2*pi*rand);
        if isfield(obj(i),'rcs')
            a = a * sqrt(max(obj(i).rcs,0));
        end

        taps(end+1)  = n; %#ok<AGROW>
        gains(end+1) = a; %#ok<AGROW>
    end

    % ---- Build FIR ----
    L = max(taps) + 1;
    h = zeros(L,1);
    for k = 1:numel(taps)
        h(taps(k)+1) = h(taps(k)+1) + gains(k);
    end
end