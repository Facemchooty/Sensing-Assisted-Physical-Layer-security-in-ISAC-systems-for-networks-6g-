function [beamAz, bestU] = selectBeamToUsers(bsPos, users, azGrid, beamwidthDeg, Gmain_dB, Gside_dB)
%SELECTBEAMTOUSERS Choose beam pointing angle that maximizes gain to one of the users
% Returns beamAz (deg) chosen from azGrid and the selected user index bestU.

    Nu = numel(users);
    bestScore = -inf;
    bestU = 1;
    beamAz = azGrid(1);

    for u = 1:Nu
        az_u = atan2d(users(u).pos(2)-bsPos(2), users(u).pos(1)-bsPos(1));

        % find best beam in codebook for this user
        scores = zeros(numel(azGrid),1);
        for k=1:numel(azGrid)
            scores(k) = sectorBeamGain(azGrid(k), az_u, beamwidthDeg, Gmain_dB, Gside_dB);
        end

        [score_u, idx] = max(scores);

        % you can bias by distance as well (optional)
        d = norm(users(u).pos - bsPos);
        score_u = score_u - 20*log10(max(d,1e-3));  % prefer closer users slightly

        if score_u > bestScore
            bestScore = score_u;
            bestU = u;
            beamAz = azGrid(idx);
        end
    end
end