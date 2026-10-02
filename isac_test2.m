% ======================================================
% Integrated Sensing and Communication (ISAC) Demo
% MVDR Beamforming + 5G-like waveform + Sensing visualization
% ======================================================

clear; close all; clc;

%% === Parameters ===
c = 3e8;                       % Speed of light
fc = 28e9;                     % 5G mmWave carrier
lambda = c/fc;
N = 16;                        % Antenna elements
d = lambda/2;
array = phased.ULA('NumElements',N,'ElementSpacing',d);

angle_user = 25;               % Communication target
angle_target = 0;              % Sensing target
angle_interf = -40;            % Interferer

SNR_dB = 20;
noiseVar = 10^(-SNR_dB/10);

%% === Simple OFDM-like waveform (toolbox-free) ===
Nfft = 512;          % FFT size
numSym = 14;         % OFDM symbols
numSC = 300;         % active subcarriers
cpLen = 64;          % cyclic prefix length

% Random QPSK data
data = pskmod(randi([0 3], numSC, numSym),4,pi/4);

% Center in frequency
ofdmSym = ifftshift([zeros((Nfft-numSC)/2, numSym); data; zeros((Nfft-numSC)/2, numSym)],1);

% OFDM modulation + cyclic prefix
txSym = ifft(ofdmSym, Nfft, 1);
cp = txSym(end-cpLen+1:end,:);
nrWaveform = reshape([cp; txSym],[],1);

fs = 122.88e6;       % typical 5G sampling rate


%% === Steering vectors ===
sv = phased.SteeringVector('SensorArray',array,...
    'PropagationSpeed',c,'IncludeElementResponse',false);

a_user   = sv(fc,angle_user);
a_target = sv(fc,angle_target);
a_int    = sv(fc,angle_interf);

%% === Covariance matrix for MVDR ===
R = (a_user*a_user') + 0.3*(a_int*a_int') + noiseVar*eye(N);
Rinv = inv(R + 1e-6*eye(N));

w_mvdr = (Rinv * a_user) / (a_user' * Rinv * a_user);
w_mvdr = w_mvdr / norm(w_mvdr);

%% === Beampattern visualization ===
angles = -90:0.5:90;
pattern = zeros(size(angles));
for k=1:length(angles)
    av = sv(fc,angles(k));
    pattern(k) = 20*log10(abs(w_mvdr' * av) + eps);
end

figure('Position',[100 100 900 400]);
subplot(1,2,1);
plot(angles,pattern,'b','LineWidth',1.4); grid on;
xlabel('Angle (deg)'); ylabel('Gain (dB)');
title('MVDR Beampattern (5G ISAC Simulation)');
xlim([-90 90]);
ylim([max(pattern)-50 max(pattern)+5]);
hold on;
plot(angle_user,interp1(angles,pattern,angle_user),'ro','MarkerFaceColor','r');
plot(angle_target,interp1(angles,pattern,angle_target),'ks','MarkerFaceColor','k');
legend('Beampattern','User','Target');

subplot(1,2,2);
theta = deg2rad(angles);
polarplot(theta, pattern - max(pattern), 'b','LineWidth',1.2);
title('Polar Beam Response');

%% === Radar-like sensing (PMCW waveform) ===
PRF = 48e3;  % Must divide fs evenly
samplesPerChip = 10;
chipWidth = samplesPerChip / fs;  % ensures integer number of samples per chip
numChips = 128;

pmcw = phased.PhaseCodedWaveform( ...
    'SampleRate', fs, ...
    'ChipWidth', chipWidth, ...
    'NumChips', numChips, ...
    'PRF', PRF, ...
    'Code', 'Barker', ...
    'NumPulses', 8);

txSig = pmcw();

% Simulate reflection from target at 50 m
tau = 2*50/c;
delaySamp = round(tau * fs);
rxSig = circshift(txSig, delaySamp) * 0.8;

% Range processing
rngProc = phased.RangeResponse( ...
    'RangeMethod', 'Matched filter', ...
    'SampleRate', fs, ...
    'PropagationSpeed', c);

[resp, rngGrid] = rngProc(txSig, rxSig);

figure('Position',[200 100 600 400]);
plot(rngGrid, 10*log10(abs(resp)/max(abs(resp))), 'LineWidth', 1.3);
xlabel('Range (m)'); ylabel('Normalized Power (dB)');
title('Range Response from PMCW Sensing (Corrected)');
grid on;



%% === Integration concept visualization ===
figure('Position',[200 100 700 500]); hold on; axis equal; grid on;
title('Integrated Sensing and Communication Visualization');
xlabel('X (m)'); ylabel('Y (m)');
xlim([-20 20]); ylim([0 25]);

Tx = [0 0];
User = [20*sind(angle_user) 20*cosd(angle_user)];
Target = [20*sind(angle_target) 20*cosd(angle_target)];
Interf = [20*sind(angle_interf) 20*cosd(angle_interf)];

plot(Tx(1),Tx(2),'kp','MarkerFaceColor','y','MarkerSize',12);
plot(User(1),User(2),'ro','MarkerFaceColor','r','MarkerSize',10);
plot(Target(1),Target(2),'ks','MarkerFaceColor','k','MarkerSize',10);
plot(Interf(1),Interf(2),'md','MarkerFaceColor','m','MarkerSize',10);
text(-1,-0.5,'Tx','FontWeight','bold');

beam = plot([Tx(1) User(1)],[Tx(2) User(2)],'g','LineWidth',2);
pause(0.5);
title('Step 1: Beamforming toward communication user');

pause(1);
delete(beam);
beam = plot([Tx(1) Target(1)],[Tx(2) Target(2)],'c--','LineWidth',2);
title('Step 2: Beam redirected for sensing (radar mode)');

pause(1);
delete(beam);
beam = plot([Tx(1) User(1)],[Tx(2) User(2)],'g','LineWidth',2);
title('Step 3: Joint ISAC mode – alternating beams for sensing & comm');
