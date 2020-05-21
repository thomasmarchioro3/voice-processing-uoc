function [c_G, G_q, eq_G] = quantize_G(G, paramfile)
% quantization
if nargin <= 1
    paramfile = 'SQ.mat';
end

    
G_range = importdata(paramfile);

dist = abs(G_range' - G);
[~, c_G] = min(dist); % map to the codeword that minimizes the error

G_q = G_range(c_G);
eq_G = abs(G - G_q);
end