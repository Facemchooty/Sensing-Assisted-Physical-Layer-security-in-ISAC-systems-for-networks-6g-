% mvdr beamforming visual
% MVDR beamforming με οπτική αναπαράσταση σε 2D & beampattern


clear; close all; clc;

%% === Array parameters ===
c = 3e8;
fc = 2.4e9;                 
lambda = c/fc;
N = 12;                     
d = lambda/2;               

array = phased.ULA('NumElements',N,'ElementSpacing',d);

%% === Geometry (angles in degrees) ===
angle_user = 25;            
angle_obstacle = 0;         
angle_interf = -40;         

SNR_dB = 20;                
INR_dB = 40;                
noiseVar = 10^(-SNR_dB/10);

sv_fun = phased.SteeringVector('SensorArray',array,'PropagationSpeed',c,'IncludeElementResponse',false);

a_user = sv_fun(fc,angle_user);
a_ob   = sv_fun(fc,angle_obstacle);
a_int  = sv_fun(fc,angle_interf);

sigma_user = 10^(SNR_dB/10);
sigma_int  = 10^((SNR_dB - 20)/10); 
sigma_ob   = 10^(INR_dB/10);          

R = sigma_user*(a_user*a_user') + ...
    sigma_int*(a_int*a_int') + ...
    sigma_ob*(a_ob*a_ob') + ...
    noiseVar*eye(N);

Rinv = inv(R + 1e-6*eye(N));
w_mvdr = (Rinv * a_user) / (a_user' * Rinv * a_user);
w_mvdr = w_mvdr / norm(w_mvdr);

%% === Array response pattern ===
angles = -90:0.5:90;
pattern_mvdr = zeros(size(angles));
for k=1:length(angles)
    av = sv_fun(fc,angles(k));
    pattern_mvdr(k) = 20*log10(abs(w_mvdr' * av) + eps);
end

%% === 1. Cartesian beam pattern plot ===
figure('Position',[100 100 900 500]);

subplot(1,2,1);
plot(angles, pattern_mvdr, 'b','LineWidth',1.4); grid on; hold on;
xlabel('Angle (deg)'); ylabel('Array response (dB)');
title('MVDR Beam Pattern (Linear View)');
xlim([-90 90]);
ylim([max(pattern_mvdr)-50 max(pattern_mvdr)+5]);

plot(angle_user, interp1(angles,pattern_mvdr,angle_user),'ro','MarkerFaceColor','r');
text(angle_user+2, interp1(angles,pattern_mvdr,angle_user),'User','Color','r');

plot(angle_obstacle, interp1(angles,pattern_mvdr,angle_obstacle),'ks','MarkerFaceColor','k');
text(angle_obstacle+2, interp1(angles,pattern_mvdr,angle_obstacle),'Obstacle','Color','k');

plot(angle_interf, interp1(angles,pattern_mvdr,angle_interf),'md','MarkerFaceColor','m');
text(angle_interf+2, interp1(angles,pattern_mvdr,angle_interf),'Interf','Color','m');

%% === 2. Polar visualization of geometry ===
subplot(1,2,2);
theta = deg2rad(angles);
r = pattern_mvdr - max(pattern_mvdr);
r = r - min(r);   % shift to positive
polarplot(theta, r, 'b', 'LineWidth',1.4); hold on;

% Mark directions
polarplot(deg2rad(angle_user), max(r),'ro','MarkerFaceColor','r');
polarplot(deg2rad(angle_obstacle), max(r)*0.95,'ks','MarkerFaceColor','k');
polarplot(deg2rad(angle_interf), max(r)*0.9,'md','MarkerFaceColor','m');
title('Spatial Visualization (Polar View)');

legend('Beam pattern','User','Obstacle','Interference','Location','southoutside');

%% === 3. 2D spatial animation (optional) ===
% Show Tx array, obstacle, and direction of main lobe

figure('Position',[200 100 800 600]);
hold on; axis equal; grid on;
xlim([-10 10]); ylim([-2 10]);
xlabel('X (m)'); ylabel('Y (m)');
title('Beamforming visualization in 2D space');

Tx = [0 0];
Rx = [8*sind(angle_user) 8*cosd(angle_user)];
Obstacle = [8*sind(angle_obstacle) 8*cosd(angle_obstacle)];
Interf = [8*sind(angle_interf) 8*cosd(angle_interf)];

plot(Tx(1),Tx(2),'kp','MarkerFaceColor','y','MarkerSize',12);
plot(Rx(1),Rx(2),'ro','MarkerFaceColor','r','MarkerSize',10);
plot(Obstacle(1),Obstacle(2),'ks','MarkerFaceColor','k','MarkerSize',10);
plot(Interf(1),Interf(2),'md','MarkerFaceColor','m','MarkerSize',10);
text(-1,-0.5,'Tx','FontWeight','bold');

% Initial beam (toward obstacle)
beam_main = plot([Tx(1) 8*sind(0)], [Tx(2) 8*cosd(0)], 'g--','LineWidth',1.5);
signalDot = plot(Tx(1),Tx(2),'ro','MarkerFaceColor','r');

pause(0.5);
title('Initial transmission: beam hits obstacle');

for k=1:40
    t = k/40;
    pos = Tx + t*(Obstacle - Tx);
    set(signalDot,'XData',pos(1),'YData',pos(2));
    pause(0.03);
end

pause(0.5);
title('Beamforming adapts: redirecting toward user');

% Beam steering animation
for step=1:50
    az = 0 + (angle_user-0)*step/50;
    delete(beam_main);
    beam_main = plot([Tx(1) 8*sind(az)], [Tx(2) 8*cosd(az)], 'g','LineWidth',2);
    pause(0.03);
end

% Move signal along new path
for k=1:40
    t = k/40;
    pos = Tx + t*(Rx - Tx);
    set(signalDot,'XData',pos(1),'YData',pos(2));
    pause(0.03);
end
title('Signal successfully reaches the user!');
