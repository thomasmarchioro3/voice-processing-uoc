clc, close, clear;

set(0,'defaultTextInterpreter','latex');

file = 'H.22.16k.wav';

extype = 2;
%{
    EXCITATON TYPES (change 'extype')
    0: normal excitation (reconstruction)
    1: whisper voice
    2: enhanced whisper voice
    3: robot voice
    4: enhanced robot voice 
%}

[sig, Fs] = audioread(file);

Horizon  = 30;   %30ms - window length
OrderLPC = 24;   %order of LPC
Buffer   = 0;    % initialization
out = zeros(size(sig)); % initialization

Horizon = Horizon*Fs/1000;
Shift   = Horizon/2;       % frame size - step size
Win     = hanning(Horizon);  % analysis window

Lsig   = length(sig);
slice  = 1:Horizon;
tosave = 1:Shift;
Nfr    = floor((Lsig-Horizon)/Shift)+1;  % number of frames

t = 0:1/Fs:length(sig)/Fs-1/Fs;

for l = 1:1:Nfr

    sigLPC = Win.*sig(slice);
    en = sum(sigLPC.^2);
    
    [r, lags] = xcorr(sigLPC);  %autocorrelation
    r(lags<0) = [];             %discarding negatives

    [a, e, k] = my_levinson(r, OrderLPC);

    G = sqrt(e);        %gain

    ex = filter(a, G, sigLPC); %inverse filter to get exitation
    ex_old = ex;
    
    switch extype
        case 1
            tname = 'whisper';
            ex = randn(Horizon, 1); 
        case 2
            tname = 'enhanced whisper';
            ex = ex + 1e-1*randn(Horizon, 1); 
        case 3
            tname = 'robot';
            ex = ones(Horizon, 1);
        case 4
            tname = 'enhanced robot';
            ex = ex + 5e-1*ones(Horizon, 1);
        otherwise
            tname = 'original';
    end
    
    ex = ex/(ex'*ex); % normalization  
    
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
    
    if l==floor(Nfr*0.33) % plots for a frame in the middle of the signal
        figure(1); 
        set(gcf,'Position', [500, 300, 420, 260]);
        plot(t(slice), ex_old); 
        hold on;
        plot(t(slice), ex);
        hold off;
        grid;
        legend('Original', 'Modified');
        xlabel('Time ($s$)');
        title(['Excitation: ', tname]);
        
        
        figure(2); 
        set(gcf,'Position', [500, 300, 420, 260]);
        plot(t(slice), sig(slice));
        hold on;
        plot(t(slice), s);
        hold off;
        grid;
        legend('Original', 'Modified');
        xlabel('Time ($s$)');
        title(['Modification of a single frame: ', tname]);
        
    end

end

figure(3);
subplot(2,1,1); 
plot(t, sig); 
legend('Original');
xlim([t(1), t(end)]);
grid;
title(['Modification of the complete signal: ', tname]);
subplot(2,1,2); plot(t, out, 'r'); 
legend('Modified'); 
xlim([t(1), t(end)]);
grid;
xlabel('Time ($s$)');


soundsc(out, Fs);