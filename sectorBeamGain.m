function GdB = sectorBeamGain(beamAz, targetAz, bwDeg, Gmain_dB, Gside_dB)
%SECTORBEAMGAIN Simple 2D sector antenna pattern
%
% beamAz     : beam pointing direction (deg)
% targetAz   : direction to user/object (deg)
% bwDeg      : beamwidth (deg)
% Gmain_dB   : mainlobe gain (dB)
% Gside_dB   : sidelobe gain (dB)

    d = wrapTo180(targetAz - beamAz);

    if abs(d) <= bwDeg/2
        GdB = Gmain_dB;
    else
        GdB = Gside_dB;
    end
end