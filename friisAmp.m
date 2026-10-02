function amp = friisAmp(d, fc)
%FRIISAMP Friis free-space amplitude (not power)
    c = 3e8;
    amp = (c/(4*pi*fc*max(d,1e-6)));
end