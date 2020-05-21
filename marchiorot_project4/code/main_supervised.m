clc, close, clear;

N_coeff = 25;
get_data('train', N_coeff);
get_data('test', N_coeff);

train = importdata('train_data.mat');
test = importdata('test_data.mat');
n_classes = length(train);

dim = size(train(1).coeff, 2);
mu = zeros(dim, n_classes);
Sigma = zeros(dim, dim, n_classes);

%apriori = ones(n_classes, 1)/n_classes;
n_coeffs = zeros(n_classes, 1);
for kk = 1:n_classes
    id19 = train(kk).coeff;
    mu(:, kk) = mean(id19);
    Sigma(:, :, kk) = cov(id19);
    n_coeffs(kk) = size(id19, 1);
end
apriori = n_coeffs/sum(n_coeffs);


gm = gmdistribution(mu', Sigma, apriori);

y_pred = predictGMM(gm, train);
y_true = find(y_pred);
accuracy = mean(y_true == y_pred);
disp(['Accuracy on training data: ' num2str(accuracy)]);

y_pred = predictGMM(gm, test);
y_true = find(y_pred);
accuracy = mean(y_true == y_pred);
disp(['Accuracy on test data: ' num2str(accuracy)]);

function [y_pred] = predictGMM(gm, data)
n_classes = gm.NumComponents;
y_pred = zeros(n_classes, 1);
for ii = 1:n_classes
    P = posterior(gm, data(ii).coeff);
    [~, labels] = max(P');
    hh = sum(labels(:) == 1:n_classes);
    [~, y_pred(ii)] = max(hh);
end
end


