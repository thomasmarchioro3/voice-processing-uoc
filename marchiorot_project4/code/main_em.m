clc, close, clear;

maxiters = 50;
N_coeff = 25;
get_data('train', N_coeff);
get_data('test', N_coeff);

train = importdata('train_data.mat');
test = importdata('test_data.mat');
n_classes = length(train);

dim = size(train(1).coeff, 2);
X = [];
for kk = 1:n_classes
    X = [X; train(kk).coeff];
end

[mu, Sigma, apriori, labels] = gmm_em(X, n_classes, maxiters);
gm = gmdistribution(mu', Sigma, apriori);

save('model_em', 'gm');