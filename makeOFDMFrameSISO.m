function tx = makeOFDMFrameSISO(Nfft,cpLen,Nsym,fftBins,M)
%MAKEOFDMFRAMESISO Generate one SISO OFDM frame (time-domain + bits)
%
% tx.time : column vector time-domain samples with CP
% tx.bits : transmit bits (column)

    bitsPerSym = log2(M);
    Nsc = numel(fftBins);

    tx.bits = randi([0 1], Nsc*Nsym*bitsPerSym, 1);

    symIdx = bi2de(reshape(tx.bits,bitsPerSym,[]).', 'left-msb');
    sym    = qammod(symIdx, M, 'UnitAveragePower', true);
    grid   = reshape(sym, Nsc, Nsym);

    X = zeros(Nfft, Nsym);
    X(fftBins,:) = grid;

    x = ifft(ifftshift(X,1), Nfft, 1);
    xcp = [x(end-cpLen+1:end,:); x];
    tx.time = xcp(:);
end