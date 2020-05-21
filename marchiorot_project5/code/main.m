clc, close, clear;
seed = 1164;
rng(seed);

filename = 'furelise-1000z';
%filename = 'marchiorot-fish';

filename_original = strcat(filename, '.wav');
filename_noisy = strcat(filename, '-noise.wav');
[s_original, Fs] = audioread(filename_original);
[s_noisy, ~] = audioread(filename_noisy);

SNR = mean(s_original.^2)./mean((s_original - s_noisy).^2);
SNR_NOISY = 20*log10(SNR);

%% UNSUPERVISED

Horizon  = 30;   %30ms - window length

%[~, sigma] = estimate_noise_gauss(s_noisy(1:1000), Fs, Horizon); 
%Sb = sigma^2*ones(floor(Horizon*Fs/1000), 1);
%Sb = estimate_noise_psd(s_noisy(1:1000), Fs, Horizon);


% SPECTRAL SUBTRACTION

lambda = 30;
%lambda = 10;

out = spectral_subtraction(s_noisy, Sb, Fs, Horizon, lambda);
outname = strcat('out/', filename, '-ss.wav');
audiowrite(outname, out, Fs);
%soundsc(out, Fs);
SNR = mean(s_original.^2)./mean((s_original - out).^2);
SNR_SS = 20*log10(SNR);

% UNSUPERVISED WIENER

lambda = 50;
out = wiener_unsupervised(s_noisy, Sb, Fs, Horizon, lambda);
%soundsc(out, Fs);
outname = strcat('out/', filename, '-wu.wav');
audiowrite(outname, out, Fs);
SNR = mean(s_original.^2)./mean((s_original - out).^2);
SNR_WU = 20*log10(SNR);


%% SUPERVISED

%SNR = mean(s_original.^2)./mean((s_original - s_noisy).^2);
%SNR_train = 20*log10(SNR);
SNR_train = 10;

% NAIVE

Hm = naive_filter(s_original, SNR_train);
Y = fft(s_noisy, length(Hm));
out = ifft(Hm.*Y, length(s_noisy));
outname = strcat('out/', filename, '-na.wav');
audiowrite(outname, out, Fs);
SNR = mean(s_original.^2)./mean((s_original - out).^2);
SNR_NA = 20*log10(SNR);


% SUPERVISED WIENER

Hs = wiener(s_original, SNR_train);
Y = fft(s_noisy, length(Hs));
out = ifft(Hs.*Y, length(s_noisy));
outname = strcat('out/', filename, '-ws.wav');
audiowrite(outname, out, Fs);
SNR = mean(s_original.^2)./mean((s_original - out).^2);
SNR_WS = 20*log10(SNR);