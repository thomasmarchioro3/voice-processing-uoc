clc, close, clear;

set(0,'defaultTextInterpreter','latex');

file = 'H.22.16k.wav';
%file = 'marchiorot_fish.wav';
%file = 'emodb_m_39_snd_norm.wav';
%file = 'arctic_bdl1_snd_norm.wav';

N = 1024;

[sig, Fs] = audioread(file);

f_shift = Fs/100;

Horizon  = 30;   %30ms - window length
OrderLPC = 14;   %order of LPC
Buffer   = 0;    % initialization
out = zeros(size(sig)); % initialization

Horizon = Horizon*Fs/1000;
Shift   = Horizon/2;       % frame size - step size
Win     = hanning(Horizon);  % analysis window

Lsig   = length(sig);
slice  = 1:Horizon;
tosave = 1:Shift;
Nfr    = floor((Lsig-Horizon)/Shift)+1;  % number of frames

% define time interval
t = 0:1/Fs:length(sig)/Fs-1/Fs;
f = 0:Fs/N:Fs - Fs/N; % range of frequencies

for l = 1:1:Nfr
    
    sigLPC = Win.*sig(slice);
    en = sum(sigLPC.^2);
    
        [r, lags] = xcorr(sigLPC);  %autocorrelation
    r(lags<0) = [];             %discarding negatives

    %[a, e, k] = levinson(r, OrderLPC); %levinson coefs.
    [a, e, k] = my_levinson(r, OrderLPC);

    G = sqrt(e);        %gain

    ex = filter(a, G, sigLPC); %inverse filter to get exitation

    % MODIFY POLES
    poles = sort(roots(a), 'ComparisonMethod', 'abs');
    poles = flip(poles); % ascending order
    
    % work on the first three most important frequencies
    %{
    p_old = poles(1:6);
    idx = find(imag(p_old)~=0);
    p_old = p_old(idx); % modify only non-real poles
    p_new = p_old;
    %}
    
    j = 1; m = 1;
    idx = [];
    while j <= 3
        while m < length(poles) && imag(poles(m)) == 0
            m = m+1;
            j = j+1;
        end
        if j > 3
            break;
        end
        idx = [idx, m];
        idx = [idx, m+1];
        m = m+2;
        j = j+1;
    end
    
    p_old = poles(idx);
    p_new = p_old;
    
    p_shift = exp(+1i*2*pi*f_shift/Fs);
    for j = 1:length(p_new)
        rule = angle(p_new(j))>0;
        %rule = (angle(p_new(j))>0 && angle(p_new(j))<pi/2) || (angle(p_new(j))>-pi && angle(p_new(j))<-pi/2);
        if rule
            p_new(j) = p_old(j)*p_shift;
        else
            p_new(j) = p_old(j)/p_shift;
        end
    end
    
    poles(idx) = p_new;
    
    % compute the coefficients from the poles
    a2 = [1, -poles(1)];
    for j = 2:1:OrderLPC

        b2 = [1, -poles(j)];
        a2 = conv(a2, b2);

    end
    a2 = real(a2);
    
    s = filter(G, a2, ex);
    
    % PLOTS
    if l==floor(Nfr*0.66) % plots for a frame in the middle of the signal
        figure(1); 
        set(gcf,'Position', [500, 100, 420, 260]);
        zplane([], p_old);
        title('Original poles');
        
        figure(2); 
        set(gcf,'Position', [800, 100, 420, 260]);
        zplane([], p_new);
        title('Modified poles');
        
        H = freqz(G, a, N);
        H2 = freqz(G, a2, N);
        
        formants = angle(p_old)/pi*Fs;
        formants = formants(formants>0);
        N_form = length(formants);
        H_form = zeros(N_form,1);
        %i_form = zeros(N_form,1);
        for i = 1:N_form
            ff = formants(i);
            [~, i_form] = min(abs(ff-f));
            formants(i) = f(i_form);
            H_form(i) = H(i_form);
        end
        
        figure(3);
        set(gcf,'Position', [500, 300, 520, 220]);
        plot(f, 10*log10(abs(H)));
        hold on;
        plot(formants, 10*log10(abs(H_form)), '^');
        hold off;
        grid;
        xlabel('$f$ (Hz)');
        title('Estimated formants (dB)');
        
        
        
        figure(4);
        set(gcf,'Position', [500, 200, 520, 420]);
        subplot(2,1,1);
        plot(f, abs(H));
        grid;
        legend('Original');
        %ylabel('Original');
        title('FFT of vocal tract filter');
        subplot(2,1,2);
        plot(f, abs(H2), 'r');
        grid;
        legend('Modified');
        %ylabel('Modified');
        
        
    end
    % synthesis
    s = filter(G, a2, ex);
    ens = sum(s.^2);     % get the short-time energy of the output
    g   = sqrt(en/ens);  % normalization factor
    s   = s*g;           % energy compensation

    s(1:Shift)  = s(1:Shift) + Buffer;    % Overlap and add
    out(tosave) = s(1:Shift);             % save the first part of the frame
    Buffer      = s(Shift+1:Horizon);     % buffer the rest of the frame

    slice  = slice + Shift;   % move the frame
    tosave = tosave + Shift;
end

S = fft(sig, N); % fft of original frame
S_tilde = fft(out, N); % fft of reconstructed frame

%{
figure(4);
subplot(2,1,1);
plot(f, 10*log10(abs(S)));
xlim([0, Fs/2]);
ylabel('$S(f)$');
subplot(2,1,2);
plot(f, 10*log10(abs(S_tilde)));
xlim([0, Fs/2]);
ylabel('$\tilde{S}(f)$');
%}

figure(5);
set(gcf,'Position', [500, 300, 520, 220]);
plot(f, 10*log10(abs(S)));
hold on;
plot(f, 10*log10(abs(S_tilde)));
hold off;
legend('Original', 'Modified');
grid;
xlim([0, Fs/2]);
title('FFT of complete signal with modified poles');


soundsc(out, Fs);