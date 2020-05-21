clc, close, clear;

N = 10000;

dim = 2;

p1 = 0.35;
p2 = 0.25;
p3 = 0.4;

X1 = 6*randn(dim,p1*N) + [12, -3]';
X2 = 4*randn(dim,p2*N) + [0, 0]';
X3 = 6*randn(dim,p3*N) + [-9, 6]';
X = [X1, X2, X3];

k = 3;

[mu, Sigma, P, labels] = gmm_em(X, k);

figure(1);
for kk = 1:k
    %color = rand(1, 3);
    X_kk = X(:, labels == kk);
    n_kk = size(X_kk, 2);
    mu_kk = mu(:, kk);
    Sigma_kk = Sigma(:, :, kk);
    p = mvnpdf(X_kk',mu_kk',Sigma_kk);
    p = p*P(kk);
    plot3(X_kk(1, :), X_kk(2, :), p, '.');
    hold on;
end
hold off;

