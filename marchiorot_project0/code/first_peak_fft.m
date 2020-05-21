function [f_peak, X_peak] = first_peak_fft(X,f)
Y = abs(X);
thresh = max(Y)*0.1;
i = 2;
while ~(Y(i)-Y(i-1)>0 && Y(i)-Y(i+1)>0 && Y(i)>thresh) 
    i = i+1;
end
f_peak = f(i);
X_peak = X(i);
end