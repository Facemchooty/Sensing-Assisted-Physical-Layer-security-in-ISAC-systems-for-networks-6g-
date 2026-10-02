% hex network beamforming


clear; close all; clc;
rng(2);

%% Παράμετροι δικτύου
numRings = 2;           % πόσοι γύροι γύρω από κέντρο (2 -> 19 cells, προσαρμόζεις)
cellRadius = 600;       % μέτρα (από κέντρο κυψέλης σε κορυφή)
numUsersPerCell = 30;   % μέγιστος χρήστες ανά κυψέλη
obsProb = 0.4;          % πιθανότητα ένας obstacle να υπάρχει σε μια κυψέλη
obsRadius = 120;        % μέγεθος εμποδίου (m)
% Array / beamforming
c = 3e8; fc = 2.4e9; lambda = c/fc;
N = 12; d = lambda/2;
array = phased.ULA('NumElements',N,'ElementSpacing',d);
sv = phased.SteeringVector('SensorArray',array,'PropagationSpeed',c,'IncludeElementResponse',false);

%% Δημιουργία hex grid κέντρων (axial coordinates)
centers = [];
for q = -numRings:numRings
    rmin = max(-numRings, -q-numRings);
    rmax = min(numRings, -q+numRings);
    for r = rmin:rmax
        x = cellRadius * (sqrt(3)*(q + r/2));
        y = cellRadius * (3/2 * r);
        centers(end+1,:) = [x,y]; %#ok<SAGROW>
    end
end
numCells = size(centers,1);

%% Για κάθε κυψέλη: users & (πιθανά) obstacle
cells = cell(numCells,1);
obstacles = nan(numCells,3); % [x,y,r] ή NaN
for i=1:numCells
    cx = centers(i,1); cy = centers(i,2);
    % τυχαίοι χρήστες εντός κυψέλης (προσέγγιση: τυχαίο πολικό εντός radius)
    nUsers = randi([floor(numUsersPerCell*0.7), numUsersPerCell]);
    th = 2*pi*rand(nUsers,1);
    rr = cellRadius*sqrt(rand(nUsers,1))*0.95; % ομοιόμορφη κατανομή εντός κύκλου προσεγγ.
    ux = cx + rr .* cos(th);
    uy = cy + rr .* sin(th);
    cells{i}.users = [ux,uy];
    cells{i}.center = [cx,cy];
    cells{i}.BS = [cx,cy]; % BS στο κέντρο
    % obstacle με κάποια πιθανότητα και τυχαία θέση εντός cell
    if rand < obsProb
        th_o = 2*pi*rand; ro = (cellRadius-obsRadius)*sqrt(rand);
        ox = cx + ro*cos(th_o); oy = cy + ro*sin(th_o);
        obstacles(i,:) = [ox,oy,obsRadius];
    end
end

%% Plot network layout
figure('Position',[100 100 900 900]);
hold on; axis equal; grid on;
title('Hex network: BS (orange), users (blue), obstacles (red triangles)');
xlabel('Meters'); ylabel('Meters');

% draw hex polygons
for i=1:numCells
    cx = centers(i,1); cy = centers(i,2);
    ang = (0:5)*pi/3 + pi/6;
    hx = cx + cellRadius*cos(ang);
    hy = cy + cellRadius*sin(ang);
    patch(hx,hy,'w','EdgeColor',[0.6 0.6 0.6]);
end

% plot users and BS
for i=1:numCells
    u = cells{i}.users;
    scatter(u(:,1), u(:,2), 8, 'b', 'filled');
    plot(cells{i}.BS(1), cells{i}.BS(2), 'p', 'MarkerSize',10, 'MarkerFaceColor',[1 .6 0], 'MarkerEdgeColor','k');
end

% plot obstacles
for i=1:numCells
    if ~isnan(obstacles(i,1))
        plot(obstacles(i,1), obstacles(i,2), 'v', 'MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',9);
        viscircles(obstacles(i,1:2), obstacles(i,3),'Color','r','LineStyle','--','LineWidth',0.5);
    end
end

%% Για κάθε BS: επιλέγουμε έναν user (τυπικά τον τυχαίο ή τον εγγύτερο) και κάνουμε MVDR
quiverScale = 0.8;
for i=1:numCells
    BSpos = cells{i}.BS;
    users = cells{i}.users;
    if isempty(users); continue; end
    % Επιλογή χρήστη: επέλεξε τυχαίο (μπορείς να αλλάξεις σε "nearest")
    idxUser = randi(size(users,1));
    userPos = users(idxUser,:);
    % Έλεγχος LOS: υπάρχει obstacle στην ίδια κυψέλη και η ευθεία περνά κοντά στο κέντρο του obstacle?
    blocked = false;
    if ~isnan(obstacles(i,1))
        obs = obstacles(i,1:2);
        r_obs = obstacles(i,3);
        % απόσταση σημείου-ευθείας
        % γραμμή από BS->user: p0 + t*(p1-p0), t in [0,1]
        p0 = BSpos; p1 = userPos;
        v = p1 - p0;
        if norm(v) > 1e-6
            t = dot((obs - p0), v) / dot(v,v);
            tclamped = max(0, min(1,t));
            closest = p0 + tclamped * v;
            dist = norm(closest - obs);
            if dist <= r_obs
                blocked = true;
            end
        end
    end
    % Υπολογισμός γωνιών (azimuth in deg) από BS
    vec = userPos - BSpos; az_user = atan2d(vec(2), vec(1));
    if blocked
        % virtual source (simple mirror): p_virtual = 2*obs - user
        o = obstacles(i,1:2);
        virt = 2*o - userPos;
        vecv = virt - BSpos; az_virt = atan2d(vecv(2), vecv(1));
        targetAz = az_virt;
    else
        targetAz = az_user;
    end
    % obstacle azimuth (if exists)
    if ~isnan(obstacles(i,1))
        o = obstacles(i,1:2);
        az_obs = atan2d((o(2)-BSpos(2)), (o(1)-BSpos(1)));
    else
        az_obs = 999;
    end
    
    % MVDR weights: steer to target, include obstacle as interferer if present
    a_target = sv(fc, targetAz);
    R = eye(N) * 1e-2; % μικρό θορυβικό baseline
    if az_obs~=999
        a_ob = sv(fc, az_obs);
        sigma_obs = 10; % ισχυρός scatterer
        R = R + sigma_obs*(a_ob*a_ob');
    end
    % small regularization
    Rinv = inv(R + 1e-6*eye(N));
    w = (Rinv*a_target) / (a_target' * Rinv * a_target);
    w = w / norm(w);
    
    % Σχεδίαση βέλους από BS προς επιδιωκόμενη γωνία
    len = cellRadius*0.6;
    dx = len * cosd(targetAz); dy = len * sind(targetAz);
    if blocked
        quiver(BSpos(1), BSpos(2), dx, dy, 'MaxHeadSize',1,'Color','m','LineWidth',1.2);
        % mark user as blocked (red)
        plot(userPos(1), userPos(2),'s','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',6);
    else
        quiver(BSpos(1), BSpos(2), dx, dy, 'MaxHeadSize',1,'Color','g','LineWidth',1.0);
        plot(userPos(1), userPos(2),'o','MarkerFaceColor','c','MarkerEdgeColor','k','MarkerSize',5);
    end
end

legend({'cell edges','','users','BS','obstacle','obs radius','beam to reflected (blocked)','blocked user','beam to user (LOS ok)'},'Location','bestoutside');
hold off;

%% Παράδειγμα: δείξε beampattern για μία επιλεγμένη κυψέλη (π.χ. την κεντρική)
% βρες κεντρική κυψέλη (με μικρότερη απόσταση από (0,0))
dists = sum(centers.^2,2);
[~, idxCenter] = min(dists);
BSpos = centers(idxCenter,:);
figure;
% επανυπολογισμός για την κεντρική κυψέλη
% διάλεξε πρώτο user της κεντρικής
userPos = cells{idxCenter}.users(1,:);
vec = userPos - BSpos; az_user = atan2d(vec(2),vec(1));
if ~isnan(obstacles(idxCenter,1))
    o = obstacles(idxCenter,1:2);
    az_obs = atan2d(o(2)-BSpos(2), o(1)-BSpos(1));
else
    az_obs = 999;
end
targetAz = az_user;
a_target = sv(fc, targetAz);
R = eye(N)*1e-2;
if az_obs~=999
    a_ob = sv(fc,az_obs);
    R = R + 10*(a_ob*a_ob');
end
Rinv = inv(R + 1e-6*eye(N));
w = (Rinv*a_target) / (a_target' * Rinv * a_target);
w = w / norm(w);

% compute pattern
angles = -180:0.5:180;
pattern_mvdr = zeros(size(angles));
for k=1:length(angles)
    av = sv(fc, angles(k));
    pattern_mvdr(k) = 20*log10(abs(w' * av) + eps);
end
plot(angles, pattern_mvdr, 'LineWidth',1.4); grid on;
xlabel('Azimuth (deg)'); ylabel('Response (dB)');
title('MVDR pattern for example (central BS)');
xlim([-180 180]);
ylim([-60 5]);
hold on;
plot(az_user, interp1(angles,pattern_mvdr,az_user),'ro','MarkerFaceColor','r');
if az_obs~=999
    plot(az_obs, interp1(angles,pattern_mvdr,az_obs),'ks','MarkerFaceColor','k');
    legend('Pattern','user az','obstacle az');
else
    legend('Pattern','user az');
end
