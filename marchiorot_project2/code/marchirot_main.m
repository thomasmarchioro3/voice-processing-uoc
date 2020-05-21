clc, close, clear;

set(0,'defaultTextInterpreter','latex');

%% file 1
file = 'arctic_bdl1_snd_norm.wav';

[X1, Fs] = audioread(file);

N = 0.030; % frame size (ms)
S = 0.015; % frame step (ms)
L = 80; % number of frequencies (sinusoids)
Ns = 32000; % number of samples of the fft

[ Y1, SNR_by_frame1 ] = SinM_marchiorot(X1, Fs, N, S, L);
SNR1 = mean(SNR_by_frame1);

soundsc(Y1, Fs);

X = fft(X1, Ns);
Y = fft(Y1, Ns);
f = Fs*(0:(Ns-1))/Ns;
plotFFT(X, Y, f, Fs);

disp('Press any key to continue... \n');

%eturn;
pause;


%% file 2


file = 'marchiorot_dont_steal.wav';

[X2, Fs] = audioread(file);

N = 0.030; % frame size (ms)
S = 0.015; % frame step (ms)
L = 80; % number of frequencies (sinusoids)

[ Y2, SNR_by_frame2 ] = SinM_marchiorot(X2, Fs, N, S, L);
SNR2 = mean(SNR_by_frame2);

soundsc(Y2, Fs);

X = fft(X2, Ns);
Y = fft(Y2, Ns);
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
legend('Original', 'Reconstructed');
subplot(2,1,2);
plot(f, angle(X));
hold on;
plot(f, angle(Y));
hold off;
grid;
xlim([0, Fs/2]);
xlabel('$f$ (Hz)');
ylabel('Phase (rad)');
legend('Original', 'Reconstructed');
%title('FFT magnitude comparison');
end