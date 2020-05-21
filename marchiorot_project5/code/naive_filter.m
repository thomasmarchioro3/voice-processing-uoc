function Hm = naive_filter(sig, SNR)

Noise = randn(size(sig));
Es = sum(sig.^2);
En = sum(Noise.^2);
sigma = sqrt(10^(-SNR/20)*Es/En);

sn = sig + sigma*Noise;
noise = sn - sig;

nn=length(sn);

X = fft(sig, nn);
N = fft(noise, nn);
Y = fft(sn, nn);

%X = Y - N = Hm.*Y = (1 - N./Y)Y

Hm = 1 - N./Y;

end