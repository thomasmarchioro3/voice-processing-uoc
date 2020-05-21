clc, close, clear;

set(0,'defaultTextInterpreter','latex');
file = 'arctic_bdl1_snd_norm.wav';

snr = 8; % snr (dB)
Ns = 32000; % number of samples of the fft

[X, Fs] = audioread(file);

N = 0.030; % frame size (ms)
S = 0.015; % frame step (ms)
L = 10; % number of frequencies (sinusoids)

% ms --> samples
N = floor(N*Fs);
N = floor(N/2)*2 + 1;
S = floor(S*Fs);

W = hanning(N);

LN = length(X);
Nfr = floor((LN - N)/S) + 1; % number of frames

Y = zeros(length(X),1);

frame = [1 N];
for i = 1:Nfr
    Xf = X(frame(1):frame(2), 1);
    sigmaN = std(Xf)*10^(-snr/20);
    R = sigmaN*randn(length(Xf), 1);
    Y(frame(1):frame(2), 1) = Xf + R;
    frame  = frame + S;
end

x = X;
y = Y;
soundsc(y, Fs);

X = fft(X, Ns);
Y = fft(Y, Ns);
f = Fs*(0:(Ns-1))/Ns;
plotFFT(X, Y, f, Fs);

function plotFFT(X, Y, f, Fs)
figure(1);
subplot(2,1,1);
plot(f, 10*log10(abs(X)));
hold on;
plot(f, 10*log10(abs(Y)));
hold off;
grid;
xlim([0, Fs/2]);
ylabel('Magnitude (dB)');
xlabel('$f$ (Hz)');
legend('Original', 'Disturbed');
subplot(2,1,2);
plot(f, angle(X));
hold on;
plot(f, angle(Y));
hold off;
grid;
xlim([0, Fs/2]);
xlabel('$f$ (Hz)');
ylabel('Phase (rad)');
legend('Original', 'Disturbed');
%title('FFT magnitude comparison');
end