function az = azimuthDeg2D(p1,p2)
%AZIMUTHDEG2D Azimuth angle in degrees from p1 to p2 in 2D
% p1, p2 are 2x1 vectors [x;y]
    v = p2 - p1;
    az = atan2d(v(2), v(1));
end