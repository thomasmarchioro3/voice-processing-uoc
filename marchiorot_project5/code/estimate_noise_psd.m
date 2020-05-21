function Sb = estimate_noise_psd(sig, Fs, Horizon)

Horizon = floor(Horizon*Fs/1000);
Shift   = floor(Horizon/2);       % frame size - step size
Win     = hanning(Horizon);  % analysis window

slice  = 1:Horizon;
Lsig   = length(sig);
Nfr    = floor((Lsig-Horizon)/Shift)+1;  % number of frames


Sb = 0;

for l = 1:1:Nfr
    sig_frame = Win.*sig(slice);
        
    Y = fft(sig_frame);
    
    Sb = Sb + Y.*conj(Y);

    slice  = slice + Shift;   % move the frame
end

Sb = Sb/Nfr;
end