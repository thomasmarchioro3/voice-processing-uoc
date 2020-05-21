function trainSQ(paramfile, b)

if nargin == 0
    paramfile = 'train_data.mat';
end
if nargin <= 1
    b = 6;
end

A = importdata(paramfile);
v_sat = A.Gmax;

L = 2^b; % number of levels;
Delta = v_sat/L; % Delta value for the range [0, v_sat]
G_range = Delta/2:Delta:v_sat-Delta/2;

save('SQ.mat', 'G_range');

end