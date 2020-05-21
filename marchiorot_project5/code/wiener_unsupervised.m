function out = wiener_unsupervised(sig, Sb, Fs, Horizon, lambda)

if nargin <= 4
    lambda = 30;
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


%Sb = sigma*ones(L, 1);

for l = 1:1:Nfr
    sig_frame = Win.*sig(slice);
    
    Y = fft(sig_frame);
    
    %L = length(Y);
    Sxhat = max(Y.*conj(Y)-lambda*Sb, 0);
    Hs = Sxhat./(Sxhat + Sb);
    
    Xhat = Hs.*Y;
        
    %s = sig_frame;
    s = real(ifft(Xhat, Horizon));
    
    s(1:Shift)  = s(1:Shift) + Buffer;    % Overlap and add
    out(tosave) = s(1:Shift);             % save the first part of the frame
    Buffer      = s(Shift+1:Horizon);     % buffer the rest of the frame
    
    slice  = slice + Shift;   % move the frame
    tosave = tosave + Shift;

end
end
