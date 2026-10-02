function bits = rxOFDMFrameSISO(y, Nfft, cpLen, Nsym, fftBins, M)
%RXOFDMFRAMESISO Simple SISO OFDM receiver (perfect timing)
%
% y      : received time-domain samples (column)
% Nfft   : FFT size
% cpLen  : cyclic prefix length
% Nsym   : number of OFDM symbols in the frame
% fftBins: active subcarrier bin indices (into 1..Nfft)
% M      : QAM order
%
% Output:
% bits: received bits (column)

    y = y(:);
    symLen = Nfft + cpLen;
    need = symLen * Nsym;

    if numel(y) < need
        error("rxOFDMFrameSISO:NotEnoughSamples", ...
              "Need %d samples but got %d. Increase frame length or check channel/filter length.", need, numel(y));
    end

    % Take exactly Nsym OFDM symbols (perfect alignment assumption)
    y = y(1:need);

    % Reshape into symbols
    Y = reshape(y, symLen, Nsym);

    % Remove CP
    Y = Y(cpLen+1:end, :);

    % FFT and shift to match ifftshift mapping at TX
    Xhat = fftshift(fft(Y, Nfft, 1), 1);

    % Extract active subcarriers
    gridHat = Xhat(fftBins, :);

    % QAM demod
    idxHat = qamdemod(gridHat(:), M, 'UnitAveragePower', true);
    bitsPerSym = log2(M);
    bits = de2bi(idxHat, bitsPerSym, 'left-msb').';
    bits = bits(:);
end