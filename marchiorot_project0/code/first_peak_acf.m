function [f_peak, acf_peak] = first_peak_acf(acf,t)

thresh = max(acf)*0.98;
i = 2;
while ~(acf(i)-acf(i-1)>0 && acf(i)-acf(i+1)>0 && acf(i)>thresh)
    if i >= length(acf)/2
        break;
    end
    i = i+1;
end

f_peak = 1/t(i);
acf_peak = acf(i);

end