clc, close, clear;

%% Parameters

frame_length_time = 30e-3;
frame_shift_time = 10e-3;

alpha = 0.6;
beta = 1.2;
gamma = 0.4;
delta = 1e-3;

%% Algorithm

%addpath('src');

[s, fs] = audioread('H.22.16k.wav');
s = s - mean(s);

[energy, ZCr] = get_metrics(s, fs, frame_length_time, frame_shift_time);

Nfr = length(energy);

% THRESHOLDS
% voiced --> high energy, low zero-crossing
% unvoiced --> low energy, high zero-crossing
% silence --> very low energy, low zero-crossing

EthreshVU = alpha*mean(energy);
EthreshUS = delta*mean(energy);
ZCRthresh = beta*mean(ZCr) - gamma*std(ZCr);

VUS = zeros(1, Nfr);

% Classification for each frame
for i = 1:1:Nfr
    if energy(i) > EthreshVU %&& ZCr(i) < ZCRthresh
        % VOICED
        VUS(i) = 1.0;
    elseif energy(i) > EthreshUS && ZCr(i) > ZCRthresh
        % UNVOICED
        VUS(i) = 0.5;
    else
        % SILENCE
        VUS(i) = 0.0;
    end
end

figure(2);
set(gcf,'Position', [500, 300, 420, 320]);
set(0, 'defaultTextInterpreter', 'latex');
plot(ZCr(VUS==1.0),energy(VUS==1.0),'d', 'DisplayName', 'Voiced');
hold on;
plot(ZCr(VUS==0.5),energy(VUS==0.5),'^', 'DisplayName', 'Unvoiced');
hold on;
plot(ZCr(VUS==0.0),energy(VUS==0.0),'o', 'DisplayName', 'Silence');
hold off;
xlim([min(ZCr),max(ZCr)]);
ylim([min(energy),max(energy)]);
grid;
legend;
title('Classification of frames on 2D plane');
xlabel('Zero Crossing');
ylabel('Energy');
%matlab2tikz('tex/classification_frames.tex');


VUSi = interp1(linspace(0,length(s), length(VUS)), VUS, 1:length(s))';

figure(3);
set(gcf,'Position', [500, 300, 820, 320]);
set(0, 'defaultTextInterpreter', 'latex');
t = 0:1/fs:length(s)/fs-1/fs;
plot(t, VUSi, 'LineWidth', 2);
hold on; 
plot(t, s/max(s), 'r'); hold off;
ylim([-1.5, 1.5]);
yticks([0,0.5,1]);
yticklabels({'Silence','Unvoiced','Voiced'})
xlabel('Time (s)');
title('Energy \& Zero-Crossings Rate-based VUS discrimination');
grid;
%matlab2tikz('tex/vus_discrimination.tex');
