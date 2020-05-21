clc, close, clear;

file = 'H.22.16k.wav';
%file = 'emodb_f_107_snd_norm.wav';
%file = 'Christine_01_neutre_snd_norm.wav';

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

    % modify the poles related to the 3 formants
    poles = sort(roots(a), 'ComparisonMethod', 'abs');
    poles = flip(poles);
    
    p_old = poles(1:6);
    p = p_old;
    d = 0; % number of discarded poles
    for i = 1:6
        if i > length(p)
            break;
        end
        if imag(p(i))
            p = p(1:end-1);
            d = d + 1;
        end
    end
    
    f_shift = Fs/4;
    poles(1:6-d) = p;
    
    a_new = [1, -poles(1)];
    for j = 2:OrderLPC           
        b = [1, -poles(j)];
        a_new = conv(a_new, b);
    end
    
    a_new = real(a_new);
    
    % synthesis
    s = filter(G, a_new, ex);
    ens = sum(s.^2);     % get the short-time energy of the output
    g   = sqrt(en/ens);  % normalization factor
    s   = s*g;           % energy compensation

    s(1:Shift)  = s(1:Shift) + Buffer;    % Overlap and add
    out(tosave) = s(1:Shift);             % save the first part of the frame
    Buffer      = s(Shift+1:Horizon);     % buffer the rest of the frame

    slice  = slice + Shift;   % move the frame
    tosave = tosave + Shift;
    
    if show_plots && l==floor(Nfr/3) % plots for a frame in the middle of the signal
        figure(1);
        zplane([], p_old); 
        title('Relevant poles');
        
        figure(2); 
        zplane([], p); 
        title('Modified relevant poles');
        
    end

end

if show_plots
    figure(3);
    subplot(2,1,1); plot(t, sig); legend('Original signal');xlim([t(1), t(end)]);
    title('Reconstruction of the complete signal');
    subplot(2,1,2); plot(t, out, 'r'); legend('Reconstructed signal'); xlim([t(1), t(end)]);
end

soundsc(out, Fs);