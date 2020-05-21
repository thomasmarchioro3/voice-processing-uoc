clc, close, clear;

set(0,'defaultTextInterpreter','latex');



file = 'arctic_bdl1_snd_norm.wav';

[X, Fs] = audioread(file);

N = 0.030; % frame size (ms)
S = 0.015; % frame step (ms)

L0 = 10;
K = 10;
deltaL = 10;
rangeL = L0:deltaL:L0+(K-1)*deltaL;

SNR_frame = zeros(1, K);
SNR = zeros(1, K);

for i = 1:K
    L = rangeL(i);
    [ Y, SNR_by_frame ] = SinM_marchiorot(X, Fs, N, S, L);
    if length(X) > length(Y)
        Y = [Y; zeros(length(X) - length(Y), 1)]; % zero padding
    elseif length(X) < length(Y)
        Y = Y(1:length(X));
    end
    SNR_frame(i) = mean(SNR_by_frame);
    SNR(i) = 20*log10(std(X)/std(X-Y));
end

figure(1);
set(gcf,'Position', [500, 300, 420, 320]);
plot(rangeL, SNR_frame, '^-');
grid;
xlabel('Number of sinusoids ($L$)');
ylabel('Average SNR for a single frame');
xlim([min(rangeL), max(rangeL)]);

figure(2);
set(gcf,'Position', [500, 300, 420, 320]);
plot(rangeL, SNR, '^-');
grid;
xlabel('Number of sinusoids ($L$)');
ylabel('Overall SNR');
xlim([min(rangeL), max(rangeL)]);


