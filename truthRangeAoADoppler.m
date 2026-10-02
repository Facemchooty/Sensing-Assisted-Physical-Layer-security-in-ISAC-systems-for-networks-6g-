function [rng, aoaDeg, fdBistatic, fdMonostatic] = truthRangeAoADoppler(bsPos, obj, lambda)
%TRUTHRANGEAOADOPPLER 2D truth range/AoA/Doppler from geometry
% bsPos: 2x1 [x;y]
% obj(i).pos: 2x1, obj(i).vel: 2x1
% lambda: wavelength
%
% Outputs (1xNo):
% rng (m), aoaDeg (deg), fdBistatic (Hz), fdMonostatic (Hz)

    No = numel(obj);
    rng = zeros(1,No);
    aoaDeg = zeros(1,No);
    fdBistatic = zeros(1,No);
    fdMonostatic = zeros(1,No);

    for i=1:No
        v = obj(i).pos - bsPos;
        d = norm(v);
        u = v / max(d,1e-12);

        rng(i) = d;
        aoaDeg(i) = atan2d(u(2), u(1));          % azimuth θ

        vRad = dot(obj(i).vel, u);               % radial velocity (m/s)
        fdBistatic(i)  = vRad / lambda;          % one-way Doppler (approx)
        fdMonostatic(i)= 2*vRad / lambda;        % two-way Doppler (radar/ISAC echo)
    end
end