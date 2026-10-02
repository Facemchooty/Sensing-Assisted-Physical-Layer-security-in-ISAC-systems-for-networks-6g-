% beamforming
clear; close all; clc;

%% --- Parameters ---
lambda = 1;                % wavelength (arbitrary units)
d = lambda/2;              % antenna spacing
N = 8;                     % number of array elements (Tx)
c = 3e8;                   % speed of light
fc = 2.4e9;
f = fc;
steerDelay = 1.2;          % seconds before beam steers away
numFrames = 150;           % total frames
pauseTime = 0.03;

% Positions
Tx = [0, 0];
Rx = [20, 0];
Obstacle = [10, 0];       % obstacle in middle
obsRadius = 1.5;

% Reflection path (mirror)
reflectPt = [10, 5];      % point above obstacle (reflection path)

%% --- Steering vector function ---
sv = @(theta) exp(1j*2*pi*(0:N-1)'*d*sind(theta)/lambda);

%% --- Figure setup ---
figure('Position',[100 100 900 600]);
axis equal; hold on; grid on;
xlim([-2 22]); ylim([-5 10]);
xlabel('X (m)'); ylabel('Y (m)');
title('Beamforming obstacle avoidance demo');

% Draw elements
plot(Tx(1),Tx(2),'kp','MarkerFaceColor','y','MarkerSize',12);
plot(Rx(1),Rx(2),'kp','MarkerFaceColor','c','MarkerSize',12);
viscircles(Obstacle,obsRadius,'Color','r','LineStyle','--');

text(Tx(1)-1,Tx(2)-1,'Transmitter');
text(Rx(1)+0.5,Rx(2),'Receiver');
text(Obstacle(1)-2,Obstacle(2)+1.2,'Obstacle');

beamLine = plot([Tx(1) Rx(1)], [Tx(2) Rx(2)], 'g-', 'LineWidth', 2);
signalDot = plot(Tx(1), Tx(2), 'ro','MarkerFaceColor','r');
reflectionLine = plot([Tx(1) reflectPt(1) Rx(1)], [Tx(2) reflectPt(2) Rx(2)], 'm--','LineWidth',1.5,'Visible','off');

%% --- Animation: initial direct beam hitting obstacle ---
fprintf('Phase 1: direct path -> blocked by obstacle\n');
for k = 1:numFrames/2
    % move "signal" along direct line
    t = k / (numFrames/2);
    pos = Tx + t * (Rx - Tx);
    set(signalDot,'XData',pos(1),'YData',pos(2));
    
    % when near obstacle, flash it
    if norm(pos - Obstacle) < obsRadius+0.3
        set(beamLine,'Color','r','LineWidth',2);
        title('Signal blocked by obstacle!');
    end
    pause(pauseTime);
end

%% --- Phase 2: adaptive beamforming (redirect path) ---
fprintf('Phase 2: beamforming redirects around obstacle\n');
set(reflectionLine,'Visible','on');
set(beamLine,'Visible','off');
title('Beamforming redirects signal around obstacle (reflection path)');
for k = 1:numFrames/2
    t = k / (numFrames/2);
    % move along Tx -> reflection point -> Rx
    if t < 0.5
        pos = Tx + 2*t * (reflectPt - Tx);
    else
        pos = reflectPt + 2*(t-0.5)*(Rx - reflectPt);
    end
    set(signalDot,'XData',pos(1),'YData',pos(2));
    pause(pauseTime);
end

title('Signal successfully reached receiver via beamforming');
fprintf('Demo complete.\n');
