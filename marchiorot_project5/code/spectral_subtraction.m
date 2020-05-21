function out = spectral_subtraction(sig, Sb, Fs, Horizon, lambda)

if nargin == 2
    Horizon = 30;
end
if nargin <= 3
    lambda = 200;
end

Buffer   = 0;    % initialization
out = zeros(size(sig)); % initialization

Horizon = floor(Horizon*Fs/1000);
Shift   = floor(Horizon/2);       % frame size - step size
Win     = hanning(Horizon);  % analysis window

Lsig   = length(sig);
slice  = 1:Horizon;
tosave = 1:Shift;
Nfr    = floor((Lsig-Horizon)/Shift)+1;  % number of frames

for l = 1:1:Nfr
    sig_frame = Win.*sig(slice);
    
    Y = fft(sig_frame);
    %L = length(Y);
    
    Xabs = max(sqrt(Y.*conj(Y)-lambda*Sb), 0);
    X = Xabs.*exp(1i*angle(Y));
    %s = sig_frame;
    s = real(ifft(X, Horizon));
    
    s(1:Shift)  = s(1:Shift) + Buffer;    % Overlap and add
    out(tosave) = s(1:Shift);             % save the first part of the frame
    Buffer      = s(Shift+1:Horizon);     % buffer the rest of the frame
    
    slice  = slice + Shift;   % move the frame
    tosave = tosave + Shift;

end

end