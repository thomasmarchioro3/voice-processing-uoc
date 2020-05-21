clc, close, clear;
% code to generate the noisy signal

% read the file: it has 1000 zero samples in the beginning
%filename = 'furelise-1000z';
filename = 'marchiorot-fish';

filename_original = strcat(filename, '.wav');
filename_noisy = strcat(filename, '-noise.wav');

[s,fs]=audioread(filename_original);

%s = s(1:10000);

% listen to the file
%soundsc(s,fs); % nice ee?

%pause;

% generate noise; gaussian - zero mean, variance is one
Noise = randn(size(s));

% set the (global) SNR in dB
SNR=10;

% find the multiplicative coeff for the noise
Es = sum(s.^2);
En = sum(Noise.^2);
a = sqrt(10^(-SNR/20)*Es/En);

sn = s+a*Noise;

% listen to the noisy file
%soundsc(sn,fs); % bad ee?

% normalize it in order to save it as wav
sn = sn/(1.1*max(abs(sn)));

SNR_true = Es/sum((s - sn).^2);
SNR_true = 20*log10(SNR_true);
% save it for your work:
audiowrite(filename_noisy,sn,fs);

% so, clean it with Spectral subtraction
