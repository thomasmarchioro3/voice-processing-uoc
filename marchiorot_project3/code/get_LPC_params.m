function [G_array, gi_array, ex_array] = get_LPC_params(file, OrderLPC, Horizon)
if nargin == 1
    OrderLPC = 5;
end
if nargin == 2
    Horizon = 30;
end

[sig, Fs] = audioread(file);
Horizon = Horizon*Fs/1000;
Shift   = Horizon/2;         % Frame size - Step size
Win     = hanning(Horizon);  % Analysis window

Lsig   = length(sig);
slice  = 1:Horizon;
Nfr    = floor((Lsig-Horizon)/Shift)+1;  % Number of frames

G_array =  zeros(1, Nfr);
gi_array =  zeros(OrderLPC, Nfr);
ex_array = zeros(Horizon, Nfr);

for l = 1:1:Nfr
    sigLPC = Win.*sig(slice);
    [r, lags] = xcorr(sigLPC);  %autocorrelation
    r(lags < 0) = [];           %discarding negatives

    [a, e, k] = my_levinson(r, OrderLPC);  % Developed function
    %k(1) = []; % remove the first 1
    G = sqrt(e);
    G_array(l) = G;
    gi_array(:, l) = log2((1-k)./(1+k));
    ex_array(:, l) = filter(a, G, sigLPC); %inverse filter to get exitation
    slice  = slice + Shift;   % move the frame
end
end