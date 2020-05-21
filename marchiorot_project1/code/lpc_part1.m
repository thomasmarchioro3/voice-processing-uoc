clc, close, clear;

set(0,'defaultTextInterpreter','latex');

file = 'H.22.16k.wav';
%file = 'emodb_f_107_snd_norm.wav';
%file = 'Christine_01_neutre_snd_norm.wav';

[sig, Fs] = audioread(file);

Horizon  = 30;   %30ms - window length
OrderLPC = 4;   %order of LPC
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

if show_plots
    t = 0:1/Fs:length(sig)/Fs-1/Fs;
end


for l = 1:1:Nfr

    sigLPC = Win.*sig(slice);
    en = sum(sigLPC.^2);
    
    [r, lags] = xcorr(sigLPC);  %autocorrelation
    r(lags<0) = [];             %discarding negatives

    %[a, e, k] = levinson(r, OrderLPC); %levinson coefs.
    [a, e, k] = my_levinson(r, OrderLPC);

    G = sqrt(e);        %gain

    ex = filter(a, G, sigLPC); %inverse filter to get exitation

    % synthesis
    s = filter(G, a, ex);
    ens = sum(s.^2);     % get the short-time energy of the output
    g   = sqrt(en/ens);  % normalization factor
    s   = s*g;           % energy compensation

    s(1:Shift)  = s(1:Shift) + Buffer;    % Overlap and add
    out(tosave) = s(1:Shift);             % save the first part of the frame
    Buffer      = s(Shift+1:Horizon);     % buffer the rest of the frame

    slice  = slice + Shift;   % move the frame
    tosave = tosave + Shift;
    
    if show_plots && l==floor(Nfr*0.33) % plots for a frame in the middle of the signal
        figure(1); 
        set(gcf,'Position', [500, 300, 420, 260]);
        plot(t(slice), ex); 
        grid;
        xlabel('Time ($s$)');
        title('Estimate frame excitation');
        
        
        figure(2); 
        set(gcf,'Position', [500, 300, 420, 260]);
        plot(t(slice), sig(slice));
        hold on;
        plot(t(slice), s);
        hold off;
        grid;
        legend('Original', 'Reconstructed');
        xlabel('Time ($s$)');
        title('Reconstruction of a single frame');
        
    end

end

figure(3);
subplot(2,1,1); 
plot(t, sig); 
legend('Original');
xlim([t(1), t(end)]);
grid;
title('Reconstruction of the complete signal');
subplot(2,1,2); plot(t, out, 'r'); 
legend('Reconstruction'); 
xlim([t(1), t(end)]);
grid;
xlabel('Time ($s$)');


soundsc(out, Fs);
