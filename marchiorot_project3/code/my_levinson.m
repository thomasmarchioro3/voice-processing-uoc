function [a, e, k] = my_levinson(r, order)

k = zeros(1, order);
q = -r(2)/r(1);
a(1) = q;

e = (1-q^2)*r(1);
k(1) = q;

for i = 2:order
    s = a*r(i:-1:2);
    q = -(r(i+1) + s)/e;
    a = [ a+q*a(i-1:-1:1), q];
    e = e*(1 - q^2);
    k(i) = q;
end

a = [1, a];
end

