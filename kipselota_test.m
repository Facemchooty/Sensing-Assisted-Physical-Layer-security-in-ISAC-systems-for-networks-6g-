% Beamforming 


clear; close all; clc;

% Παράμετροι array
c = 3e8;
fc = 2.4e9;                 % συχνότητα (Hz)
lambda = c/fc;
N = 12;                     % αριθμός στοιχείων (ULA)
d = lambda/2;               % διάστημα στοιχείων

% Δημιουργία ULA
array = phased.ULA('NumElements',N,'ElementSpacing',d);

% Γωνίες (σε μοίρες)
angle_user = 25;            % γωνία χρήστη (σε μοίρες, azimuth)
angle_obstacle = 0;         % θέση εμποδίου (π.χ. μπροστά στο 0°)
angle_interf = -40;         % πιθανή άλλη παρεμβολή

% Σήματα / ισχείς (για covariance simulation)
SNR_dB = 20;                % SNR του χρήστη
INR_dB = 40;                % ισχύς "αντικειμένου/παρεμβολής" (αν το εμπόδιο δρα σαν ισχυρός scatterer)
noiseVar = 10^(-SNR_dB/10);

% Steering vectors
sv_user = phased.SteeringVector('SensorArray',array,'PropagationSpeed',c,'IncludeElementResponse',false);
a_user = sv_user(fc,angle_user);

sv_ob = phased.SteeringVector('SensorArray',array,'PropagationSpeed',c,'IncludeElementResponse',false);
a_ob = sv_ob(fc,angle_obstacle);

sv_int = phased.SteeringVector('SensorArray',array,'PropagationSpeed',c,'IncludeElementResponse',false);
a_int = sv_int(fc,angle_interf);

% Σχηματισμός τεχνητού covariance matrix:
% μοντέλο: R = sigma_user * a_user*a_user' + sigma_int * a_int*a_int' + sigma_obstacle * a_ob*a_ob' + noise*I
sigma_user = 10^(SNR_dB/10);
sigma_int  = 10^( (SNR_dB - 20)/10 ); % μικρότερη παρεμβολή (προσαρμόστε)
sigma_ob   = 10^(INR_dB/10);          % ισχυρός scatterer/εμπόδιο που "απορροφά/σπρώχνει" ενέργεια

R = sigma_user*(a_user*a_user') + sigma_int*(a_int*a_int') + sigma_ob*(a_ob*a_ob') + noiseVar*eye(N);

% Υπολογισμός MVDR weights (κλασικός closed-form)
Rinv = inv(R + 1e-6*eye(N));  % μικρή regularization
w_mvdr = (Rinv * a_user) / (a_user' * Rinv * a_user);  % MVDR / Capon

% Κανονικοποίηση (προαιρετικη)
w_mvdr = w_mvdr / norm(w_mvdr);

% Αναπαράσταση beampattern
angles = -90:0.5:90;
bv = phased.ArrayResponse('SensorArray',array,'PropagationSpeed',c);
pattern_mvdr = zeros(size(angles));
for k=1:length(angles)
    av = sv_user(fc,angles(k));
    pattern_mvdr(k) = 20*log10( abs(w_mvdr' * av) + eps );
end

% Plot beampattern
figure;
plot(angles, pattern_mvdr, 'LineWidth',1.4);
grid on;
xlabel('Angle (deg)');
ylabel('Array response (dB)');
title('MVDR beamformer pattern');
xlim([-90 90]);
ylim([max(pattern_mvdr)-50 max(pattern_mvdr)+5]);
hold on;
% επισημάνσεις
plot(angle_user, interp1(angles,pattern_mvdr,angle_user), 'ro','MarkerFaceColor','r');
text(angle_user+2, interp1(angles,pattern_mvdr,angle_user), 'User','Color','r');
plot(angle_obstacle, interp1(angles,pattern_mvdr,angle_obstacle), 'ks','MarkerFaceColor','k');
text(angle_obstacle+2, interp1(angles,pattern_mvdr,angle_obstacle), 'Obstacle','Color','k');
plot(angle_interf, interp1(angles,pattern_mvdr,angle_interf), 'md','MarkerFaceColor','m');
text(angle_interf+2, interp1(angles,pattern_mvdr,angle_interf), 'Interf','Color','m');

% Εμφάνιση weights (φάση & μέτρο)
figure;
subplot(2,1,1);
stem(0:N-1, abs(w_mvdr),'filled');
xlabel('Element index'); ylabel('Magnitude');
title('MVDR weights magnitude');

subplot(2,1,2);
stem(0:N-1, angle(w_mvdr)*180/pi,'filled');
xlabel('Element index'); ylabel('Phase (deg)');
title('MVDR weights phase');

% Προσομοίωση σήματος (προαιρετικά): στέλνουμε a_user και βλέπουμε έξοδο
sig = 1; % μονότονο σύμβολο
x = sqrt(sigma_user)*a_user*sig + sqrt(sigma_int)*a_int*0.1 + sqrt(sigma_ob)*a_ob*0.01 + sqrt(noiseVar)*(randn(N,1)+1j*randn(N,1))/sqrt(2);
y = w_mvdr' * x;
fprintf('Output SNR approx (linear): %.2f\n', (abs(w_mvdr'* (sqrt(sigma_user)*a_user))^2) / (abs(w_mvdr'*(sqrt(sigma_int)*a_int + sqrt(sigma_ob)*a_ob))^2 + noiseVar*norm(w_mvdr)^2) );
