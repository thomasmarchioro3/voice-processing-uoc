function [c_g, g_q, eq_g] = quantize_gi(g_in, paramfile)
if nargin == 1
    paramfile = 'VQ.mat';
end
centroids = importdata(paramfile);
%centroids = centroids';
if size(g_in, 2) ~= size(centroids, 2)
    g_in = g_in';
end
g_q = zeros(size(g_in));
c_g = zeros(size(g_in, 1), 1);
eq_g = zeros(size(g_in, 1), 1);

%disp(size(g_in))


for ii = 1:size(g_in, 1)
    v = g_in(ii, :);
    dist = sum((abs(v - centroids)), 2);
    [eq, idx] = min(dist);
    g_q(ii, :) = centroids(idx, :);
    c_g(ii) = idx;
    eq_g(ii) = eq;
end