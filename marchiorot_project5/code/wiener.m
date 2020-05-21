function Hs = wiener(sig, SNR)

Noise = randn(size(sig));
Es = sum(sig.^2);
En = sum(Noise.^2);
sigma = sqrt(10^(-SNR/20)*Es/En);

sn = sig + sigma*Noise;
noise = sn - sig;

%nn=length(sig);
nn = length(sig);

X = fft(sig, nn);
N = fft(noise, nn);
Y = fft(sn, nn);

Sx = X.*conj(X);
Sn = N.*conj(N);

Hs = Sx./(Sx + Sn);

end
