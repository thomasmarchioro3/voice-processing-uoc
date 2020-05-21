clc, close, clear;

%[s, fs] = audioread('H.22.16k.wav');
%[s, fs] = audioread('aa-pout-105Hz-8kHz.wav');
%[s, fs] = audioread('aa-pout-115Hz-8kHz.wav');
%[s, fs] = audioread('ae-pout-230Hz-8kHz.wav');
[s, fs] = audioread('ae-pout-255Hz-8kHz.wav');
%[s, fs] = audioread('ae-pout-315Hz-8kHz.wav');
%s = s(100:1100);

t = 0:1/fs:length(s)/fs-1/fs;
%plot(t,s);

%% Params

N = 4096;

frame_length_time = 30e-3;
frame_shift_time = 10e-3;


%% Algorithm

D = length(s);
L = frame_length_time*fs;

U = frame_shift_time*fs;

win = hamming(L);

% Number of frames
Nfr = ceil((D-L)/U);

f_peak1 = zeros(1,Nfr);
f_peak2 = zeros(1,Nfr);

%t_frames = frame_length_time/2:frame_shift_time:frame_length_time/2+(Nfr-1)*frame_shift_time;
t1=frame_length_time/2;
t2=frame_length_time/2+Nfr*frame_shift_time-frame_shift_time;
t_frames = t1:frame_shift_time:t2;

for i = 1:1:Nfr
    frame = s((i-1)*U+1:(i-1)*U+L).*win; % windowed frame
    [f_peak1(i), f_peak2(i)] = estimate_pitch(frame,fs);    
end

pitch1 = interp1(t_frames, f_peak1, t, 'spline');
pitch2 = interp1(t_frames, f_peak2, t, 'spline');

avg_pitch = mean(f_peak1(f_peak1 <= 500 & f_peak1 >= 70));
%avg_pitch = mean(f_peak2(f_peak2 <= 500 & f_peak2 >= 70));

agesex = 'child';

if avg_pitch < 275
    % ADULT
    if avg_pitch < 155
        agesex = 'adult male';
    else
        agesex = 'adult female';
    end
end

output = ['Estimate age and sex: ' agesex];
disp(output);

figure;
subplot(211);
plot(t,pitch1/1000, 'DisplayName', 'FFT');
hold on;
plot(t,pitch2/1000, 'DisplayName', 'ACF');
hold off;
ylim([0,fs/2000]);
legend;
ylabel('Frequency (kHz)');
title(output);
subplot(212); 
spectrogram(s, L, U, N, fs, 'yaxis'); 
xlabel('Time (s)');



%% helper functions

function [f_peak1, f_peak2] = estimate_pitch(x,fs)
N = 2048;

X = fft(x, N);
acf_x = xcorr(x);

f = 1:fs/N:fs;
t = -length(x)/fs+1/fs:1/fs:length(x)/fs-1/fs;

% keep only the frequencies below fs/2
X = X(f<fs/2);
f = f(f<fs/2);

[f_peak1, ~] = first_peak_fft(X,f);

% keep only autocorrelation with respect to positive shifts (k>0)
acf_x = acf_x(t>0);
t = t(t>0);

[f_peak2, ~] = first_peak_acf(acf_x,t);
end