clc, close, clear;

%addpath('src');

%[s, fs] = audioread('H.22.16k.wav');
%[s, fs] = audioread('ae-pout-230Hz-8kHz.wav');
%[s, fs] = audioread('ae-pout-255Hz-8kHz.wav');
[s, fs] = audioread('ae-pout-315Hz-8kHz.wav');

N = 4096;

frame_length_time = 30e-3;
frame_shift_time = 10e-3;

D = length(s);
L = frame_length_time*fs;

U = frame_shift_time*fs;

win = hamming(L);

% Number of frames
Nfr = ceil((D-L)/U);

i = 13;
frame = s((i-1)*U+1:(i-1)*U+L).*win; % windowed frame

x = frame;
X = fft(frame, N);
f = 1:fs/N:fs;
acf_x = xcorr(frame);
%t = linspace(-length(x)/fs, length(x)/fs, length(acf_x));
t = -length(x)/fs+1/fs:1/fs:length(x)/fs-1/fs;

% keep only the frequencies below fs/2
X = X(f<fs/2);
f = f(f<fs/2);

[f_peak, X_peak] = first_peak_fft(X,f);
t_peak = 1/f_peak;

figure(1);
set(0, 'defaultTextInterpreter', 'latex');
plot(f,abs(X));
hold on;
plot(f_peak, abs(X_peak), '^');
hold off;
ylabel('Magnitude'); xlabel('$f$ (Hz)') 
ax = gca; ax.GridLineStyle = ':'; 
grid;
title('FFT of frame (half band)');
%matlab2tikz('tex/fft.tex');

% keep only autocorrelation with respect to positive shifts (k>0)
acf_x = acf_x(t>0);
t = t(t>0);

[f_peak2, acf_peak] = first_peak_acf(acf_x,t);
t_peak2 = 1/f_peak2;

figure(2);
set(0, 'defaultTextInterpreter', 'latex');
plot(t, acf_x);
hold on;
plot(t_peak2, acf_peak, '^');
hold off;
ylabel('Autocorrelation'); xlabel('$t$ (s)') 
ax = gca; ax.GridLineStyle = ':'; 
grid;
title('Autocorrelation of frame (positive shifts)');
%matlab2tikz('tex/autocorr.tex');


