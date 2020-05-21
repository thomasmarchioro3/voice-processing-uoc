clc, close, clear;

[sig, Fs] = audioread('marchiorot_fish.wav');
Gq_codebook = importdata('SQ.mat');
gq_codebook = importdata('VQ.mat'); gq_codebook = gq_codebook';

Horizon  = 30;   %30ms - window length
OrderLPC = 6;   %order of LPC
Buffer   = 0;    % initialization
out = zeros(size(sig)); % initialization

Horizon = Horizon*Fs/1000;
Shift   = Horizon/2;       % frame size - step size
Win     = hanning(Horizon);  % analysis window

Lsig   = length(sig);
slice  = 1:Horizon;
tosave = 1:Shift;
Nfr    = floor((Lsig-Horizon)/Shift)+1;  % number of frames

cG_array =  zeros(1, Nfr);
cg_array =  zeros(1, Nfr);
ex_array = zeros(Horizon, Nfr);

for l = 1:1:Nfr

    sigLPC = Win.*sig(slice);
    en = sum(sigLPC.^2);
    
    [r, lags] = xcorr(sigLPC);  %autocorrelation
    r(lags<0) = [];             %discarding negatives

    [a, e, k] = my_levinson(r, OrderLPC);
    g = log2((1-k)./(1+k));
    G = sqrt(e);        %gain
    ex = filter(a, G, sigLPC); %inverse filter to get exitation
    
    % quantization
    [c_G, ~, ~] = quantize_G(G, 'SQ.mat');
    [c_g, ~, ~] = quantize_gi(g, 'VQ.mat');
    
    
    cG_array(l) = c_G;
    cg_array(l) = c_g;
    ex_array(:, l) = ex;
    
    %% decoder
    
    cG = cG_array(l);
    cg = cg_array(l);
    ex_l = ex_array(:, l); 
    G_q = Gq_codebook(cG);
    g_q = gq_codebook(:, cg);
    
    k_q = ((1 - 2.^g_q)./(2.^g_q + 1))';
    a_q = latc2tf(k_q);

    % synthesis
    so = filter(G_q, a_q, ex_l);
    %ens = sum(s.^2);     % get the short-time energy of the output
    %gg   = sqrt(en/ens);  % normalization factor
    %s   = s*gg;           % energy compensation

    so(1:Shift)  = so(1:Shift) + Buffer;    % Overlap and add
    out(tosave) = so(1:Shift);             % save the first part of the frame
    Buffer      = so(Shift+1:Horizon);     % buffer the rest of the frame

    slice  = slice + Shift;   % move the frame
    tosave = tosave + Shift;

end

soundsc(out, Fs);
audiowrite('marchiorot_output.wav', out, Fs);