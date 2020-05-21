clc, close, clear;
set(0,'defaultTextInterpreter','latex');
[s, Fs] = audioread('marchiorot_fish.wav');
out = audioread('marchiorot_output.wav');

N = 10000;
t = 0:1/Fs:length(s)/Fs-1/Fs;
f = 0:Fs/N:Fs - Fs/N; % range of frequencies

% time domain
figure(1);
subplot(2,1,1);
plot(t, s, 'DisplayName', 'Original');
xlabel('Time ($s$)');
legend;
subplot(2,1,2);
plot(t, out, 'r', 'DisplayName', 'Reconstruction');
xlabel('Time ($s$)');
legend;

% frequency domain
S = fft(s, N);
O = fft(out, N);
figure(2);
subplot(2,1,1);
plot(f, 10*log10(abs(S)), 'DisplayName', 'Original');
hold on;
plot(f, 10*log10(abs(O)), 'DisplayName', 'Reconstruction');
hold off;
xlabel('$f$ (Hz)');
title('Magnitude (dB)');
legend;
subplot(2,1,2);
plot(f, angle(S), 'DisplayName', 'Original');
hold on;
plot(f, angle(O), 'DisplayName', 'Reconstruction');
hold off;
xlabel('$f$ (Hz)');
title('Phase');
legend;



