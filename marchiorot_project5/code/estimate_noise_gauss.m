function [mu, sigma] = estimate_noise_gauss(sig, Fs, Horizon)

Horizon = floor(Horizon*Fs/1000);
Shift   = floor(Horizon/2);       % frame size - step size

slice  = 1:Horizon;
Lsig   = length(sig);
Nfr    = floor((Lsig-Horizon)/Shift)+1;  % number of frames

mu = 0;
for l = 1:Nfr
    sig_frame = sig(slice);
    mu = mu + sig_frame;
    slice  = slice + Shift;   % move the frame
end
mu = mu/Nfr; % unbiased estimator of the mean

slice  = 1:Horizon;
variance = 0;
for l = 1:Nfr
    sig_frame = sig(slice);
    variance = variance + (sig_frame-mu).^2;
    slice  = slice + Shift;   % move the frame
end
sigma = sqrt(variance/Nfr); % biased estimator of stddev

mu = mean(mu);
sigma = mean(sigma);

end