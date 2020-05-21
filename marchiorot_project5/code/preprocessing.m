clc, close, clear;

filename = 'marchiorot-fish';

filename_original = strcat(filename, '.wav');
%filename_noisy = strcat(filename, '-noise.wav');

[s,fs]=audioread(filename_original);

s = s - mean(s);

s_cat = [];
while length(s_cat) < 300000
    s_cat = [s_cat; s];
end
s = s_cat;

if sum(s(1:1000).^2) > 0
    s = [zeros(1000, 1); s];
end

s = s/(1.1*max(abs(s)));

audiowrite(filename_original, s, fs);