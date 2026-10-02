function tx = makeOFDMFrame(Nfft, cpLen, Nsym, fftBins, M)
%MAKEOFDMFRAME Generate one OFDM frame (time-domain + active subcarrier grid + bits)
%
% tx.time : column vector time-domain samples with CP
% tx.grid : Nsc x Nsym QAM symbols on active subcarriers
% tx.bits : transmit bits (column)

    bitsPerSym = log2(M);
    Nsc = numel(fftBins);

    % Random bits
    tx.bits = randi([0 1], Nsc*Nsym*bitsPerSym, 1);

    % QAM modulation (requires Communications Toolbox)
    symIdx = bi2de(reshape(tx.bits, bitsPerSym, []).', 'left-msb');
    sym    = qammod(symIdx, M, 'UnitAveragePower', true);

    % Grid on active subcarriers
    tx.grid = reshape(sym, Nsc, Nsym);

    % Map to full FFT bins
    X = zeros(Nfft, Nsym);
    X(fftBins,:) = tx.grid;

    % IFFT (ifftshift for centered mapping)
    x = ifft(ifftshift(X,1), Nfft, 1);

    % Add CP and serialize
    xcp = [x(end-cpLen+1:end,:); x];
    tx.time = xcp(:);
end