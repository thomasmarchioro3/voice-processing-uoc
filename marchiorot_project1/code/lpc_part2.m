clc, close, clear;

set(0,'defaultTextInterpreter','latex'); % latex font in plots

%file = 'emodb_f_107_snd_norm.wav';
file = 'H.22.16k.wav';

[sig, Fs] = audioread(file);

N = 1024; % number of samples for the FFT

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
        % frequency analysis        
        X = fft(sig(slice), N); % fft of original frame
        X_tilde = fft(s, N); % fft of reconstructed frame
        H = freqz(G, a, N); % frequency response of the filter
        Ug = fft(ex, N); % fft of the excitation
        
        % display frequencies in range [0, fmax] for better resolution
        fmax = Fs/2; 
        %Nmax = floor(N*fmax/Fs);
        
        
        figure(1);
        subplot(2,1,1);
        plot(f, abs(X));
        grid;
        xlim([0, fmax]);
        ylabel('$X(f)$');
        legend('Original');
        title('FFT of a single frame');
        subplot(2,1,2);
        plot(f, abs(X_tilde), 'r');
        grid;
        xlim([0, fmax]);
        xlabel('$f$ (Hz)');
        ylabel('$\tilde{X}(f)$');
        legend('Reconstruction');
        
        
        figure(2);        
        subplot(2,1,1);
        plot(f, abs(H));
        grid;
        %xlim([0, fmax]);
        ylabel('$H(f)$');
        legend('Vocal tract filter');
        title('FFT of a single frame');
        subplot(2,1,2);
        plot(f, abs(Ug));
        grid;
        %xlim([0, fmax]);
        xlabel('$f$ (Hz)');
        ylabel('$U_{g}(f)$');
        legend('Excitation');
        
    end

end

S = fft(sig, N); % fft of original frame
S_tilde = fft(out, N); % fft of reconstructed frame

%{
figure(3);
subplot(2,1,1);
plot(f, abs(S));
xlim([0, fmax]);
ylabel('$S(f)$');
subplot(2,1,2);
plot(f, abs(S_tilde));
xlim([0, fmax]);
ylabel('$\tilde{S}(f)$');
%}


figure(3);
set(gcf,'Position', [500, 300, 520, 220]);
plot(f, 10*log10(abs(S)));
hold on;
plot(f, 10*log10(abs(S_tilde)));
hold off;
legend('Original', 'Reconstruction');
ylabel('FFT Magnitude (dB)')
xlabel('$f$ (Hz)');
xlim([0, fmax]);
%title('Reconstruction of the complete signal');
grid;

soundsc(out, Fs);