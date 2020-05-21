clc, close, clear;

maxiters = 50;
N_coeff = 25;
get_data('train', N_coeff);
get_data('test', N_coeff);

train = importdata('train_data.mat');
test = importdata('test_data.mat');
n_classes = length(train);

dim = size(train(1).coeff, 2);
mu = zeros(dim, n_classes);
Sigma = zeros(dim, dim, n_classes);

X = [];
for kk = 1:n_classes
    X = [X; train(kk).coeff];
end

labels = kmeans(X, n_classes);
apriori = ones(n_classes, 1)/n_classes;

for kk = 1:n_classes
    X_kk = X(labels == kk, :);
    mu(:, kk) = mean(X_kk);
    Sigma(:, :, kk) = cov(X_kk);
end

gm = gmdistribution(mu', Sigma, apriori);

save('model_kmeans', 'gm');