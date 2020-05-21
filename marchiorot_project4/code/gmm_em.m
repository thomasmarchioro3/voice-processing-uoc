function [mu, Sigma, Pqk, labels] = gmm_em(X, k, maxiters)
if nargin <= 2
    maxiters = 20;
end

if size(X,1) > size(X,2)
    X = X';
end
dim = size(X,1);
n = size(X,2);
labels = ceil(k*rand(1,n)); % initialize labels at random
Pqk = mean(labels(:) == 1:k);
mu = rand(dim, k);
Sigma = rand(dim, dim, k);
%return
for kk = 1:k
    X_kk = X(:, labels == kk);
    mu(:, kk) = mean(X_kk, 2);
    Sigma(:, :, kk) = cov(X_kk');
end
iters = 0;
condition = false;
while ~condition
    labels_old = labels;
    log_lh = zeros(k, n);
    for kk = 1:k
        mu_kk = mu(:, kk);
        Sigma_kk = Sigma(:, :, kk);
        p_x_given_qk = mvnpdf(X',mu_kk',Sigma_kk); % p(X|qk, model)
        p_qk = Pqk(kk); % p(qk|model)
        log_lh(kk, :) = log(p_x_given_qk) + log(p_qk); % compute log-likelihood approximation
    end
    [~, labels] = max(log_lh);
    
    % update the parameters
    Pqk = mean(labels(:) == 1:k);
    for kk = 1:k
        X_kk = X(:, labels == kk);
        mu(:, kk) = mean(X_kk, 2);
        Sigma(:, :, kk) = cov(X_kk');
    end
    iters = iters + 1;
    condition = sum(labels == labels_old) == n; % if no further changes occur, stop
    condition = condition || iters > maxiters;
end
end