clc, close, clear;

%file = 'marchiorot_fish.wav';
file = 'H.22.16k.wav';
[sig, Fs] = audioread(file);


%file = 'H.22.16k.wav';
file = 'marchiorot_fish.wav';
%file = 'lion_1522.wav';
[sig2, Fs2] = audioread(file);

sig = sig(1:min(length(sig), length(sig2)));
sig2 = sig2(1:length(sig));

Horizon  = 30;   %30ms - window length
OrderLPC = 10;   %order of LPC
Buffer   = 0;    % initialization
out = zeros(size(sig)); % initialization

Horizon = Horizon*Fs/1000;
Shift   = Horizon/2;       % frame size - step size
Win     = hanning(Horizon);  % analysis window

Lsig   = length(sig);
slice  = 1:Horizon;
tosave = 1:Shift;
Nfr    = floor((Lsig-Horizon)/Shift)+1;  % number of frames

show_plots = true;


t = 0:1/Fs:length(sig)/Fs-1/Fs;


for l = 1:1:Nfr

    sigLPC = Win.*sig(slice);
    en = sum(sigLPC.^2);
    
    sigLPC2 = Win.*sig2(slice);    
    
    [G1, a1, ex1] = my_lpc(sigLPC, OrderLPC);
    [G2, a2, ex2] = my_lpc(sigLPC2, OrderLPC);
    
    % synthesis
    s = filter(G1, a1, ex2); % filter from sig1, excitation from sig2
    ens = sum(s.^2);     % get the short-time energy of the output
    g   = sqrt(en/ens);  % normalization factor
    s   = s*g;           % energy compensation

    s(1:Shift)  = s(1:Shift) + Buffer;    % Overlap and add
    out(tosave) = s(1:Shift);             % save the first part of the frame
    Buffer      = s(Shift+1:Horizon);     % buffer the rest of the frame

    slice  = slice + Shift;   % move the frame
    tosave = tosave + Shift;
    
    if l==floor(Nfr*0.4) % plots for a frame in the middle of the signal
        
        figure(1); 
        set(gcf,'Position', [500, 300, 420, 260]);
        plot(t(slice), sig(slice));
        hold on;
        plot(t(slice), sig2(slice));
        hold off;
        grid;
        xlabel('Time ($s$)');
        legend('Signal 1', 'Signal 2');
        title('Frames');
        
        figure(2); 
        subplot(2,1,1);
        plot(t(slice), ex1); 
        grid;
        legend('Signal 1');
        title('Estimate frame excitations');
        subplot(2,1,2);
        plot(t(slice), ex2, 'r');
        grid;
        legend('Signal 2');
        
        
        
        figure(3); 
        set(gcf,'Position', [500, 300, 420, 260]);
        plot(t(slice), sig(slice));
        hold on;
        plot(t(slice), s);
        hold off;
        grid;
        legend('Original', 'Reconstructed');
        title('Reconstruction of a single frame');
    end

end


figure(4);
subplot(2,1,1); 
plot(t, sig); 
grid;
legend('Original');
xlim([t(1), t(end)]);
title('Complete signal with changed excitation');
subplot(2,1,2); 
plot(t, out, 'r'); 
grid;
legend('Reconstruction'); 
xlim([t(1), t(end)]);


soundsc(out, Fs);

function [G, a, ex] = my_lpc(sigLPC, OrderLPC)
    [r, lags] = xcorr(sigLPC);  %autocorrelation
    r(lags<0) = [];             %discarding negatives

    %[a, e, k] = levinson(r, OrderLPC); %levinson coefs.
    [a, e, k] = my_levinson(r, OrderLPC);
    
    if ~isempty(k(abs(k)>1))
        disp('Error: unstable system.')
    end

    G = sqrt(e);        %gain

    ex = filter(a, G, sigLPC); %inverse filter to get exitation
end