clc, close, clear;
gm = importdata('model_kmeans.mat');
%gm = importdata('model_em.mat');

train = importdata('train_data.mat');
test = importdata('test_data.mat');
n_classes = length(train);

y_pred = predictGMM(gm, train);
accuracy = length(unique(y_pred))/n_classes;
disp(['Accuracy on training data: ' num2str(accuracy)]);

y_pred = predictGMM(gm, test);
accuracy = length(unique(y_pred))/n_classes;
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