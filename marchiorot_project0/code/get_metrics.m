function [energy,ZCr]=get_metrics(s, fs, frame_length_time, frame_shift_time)

D = length(s);
% 30 ms per frame, since the sampling frequency is 16 kHz, are
% L=30e-3*16e3=480 samples per frame
L = frame_length_time*fs;

% The same reasoning holds for the frame shift
U = frame_shift_time*fs;

% Window type (Hamming)
win = hamming(L);

% Number of frames
Nfr = ceil((D-L)/U); % subtract L/U frames to avoid exceeding the array
                   % alternatively one could perform zero padding

% Memory allocation (for speed)
energy = zeros(1, Nfr+1);
ZCr = zeros(1, Nfr+1);

% Loop which calculates the speech features
for i = 1:1:Nfr
    frame = s((i-1)*U+1:(i-1)*U+L).*win; % windowed frame
    energy(i) = sum(frame.^2)/L; % calculate energy
    ZCr(i) = 0.5*sum(abs(sign(frame(2:end)) - sign(frame(1:end-1))));
end

end
