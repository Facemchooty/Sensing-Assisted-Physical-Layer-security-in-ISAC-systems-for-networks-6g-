function [rngWall, hitPoint] = rangeToRoomWalls(bsPos, beamAzDeg, room)
%RANGETOROOMWALLS Distance from BS to nearest room boundary along ray direction
% Room boundaries: x=0, x=L, y=0, y=W

    x0 = bsPos(1); y0 = bsPos(2);
    th = deg2rad(beamAzDeg);
    dx = cos(th); dy = sin(th);

    % avoid divide-by-zero
    epsv = 1e-12;
    if abs(dx) < epsv, dx = sign(dx+epsv)*epsv; end
    if abs(dy) < epsv, dy = sign(dy+epsv)*epsv; end

    candT = [];

    % Intersect with x=0 and x=L
    t1 = (0 - x0)/dx;         y1 = y0 + t1*dy;
    if t1 > 0 && y1 >= 0 && y1 <= room.W, candT(end+1) = t1; end %#ok<AGROW>

    t2 = (room.L - x0)/dx;    y2 = y0 + t2*dy;
    if t2 > 0 && y2 >= 0 && y2 <= room.W, candT(end+1) = t2; end %#ok<AGROW>

    % Intersect with y=0 and y=W
    t3 = (0 - y0)/dy;         x3 = x0 + t3*dx;
    if t3 > 0 && x3 >= 0 && x3 <= room.L, candT(end+1) = t3; end %#ok<AGROW>

    t4 = (room.W - y0)/dy;    x4 = x0 + t4*dx;
    if t4 > 0 && x4 >= 0 && x4 <= room.L, candT(end+1) = t4; end %#ok<AGROW>

    if isempty(candT)
        rngWall = NaN;
        hitPoint = [NaN; NaN];
        return;
    end

    tmin = min(candT);
    hitPoint = [x0 + tmin*dx; y0 + tmin*dy];
    rngWall = tmin;  % since direction vector is unit length
end