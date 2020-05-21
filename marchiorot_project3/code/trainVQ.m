function trainVQ(paramfile, b, N_train)

if nargin == 0
    paramfile = 'train_data.mat';
end
if nargin <= 1
    b = 13;
end

A = importdata(paramfile);
g_train = A.g_array;
if size(g_train, 2) > size(g_train, 1)
    g_train = g_train';
end

N = size(g_train, 1); % number of training samples
if nargin <= 2
    N_train = N;
end

N_train = min(N_train, N);
g_train = g_train(randperm(N, N_train), :);
L = 2^b; % number of clusters
%num_clusters = 1;

num_clusters = 1;
labels = ones(N_train, 1);
centroids = zeros(L, size(g_train, 2));
while num_clusters < L
    for jj = 1:num_clusters
        g_jj = g_train(labels == jj, :);
        if size(g_jj, 1) > 1
            [idx, c_jj] = kmeans(g_jj, 2);
            num_clusters = num_clusters + 1;
            idx(idx == 1) = jj; 
            idx(idx == 2) = num_clusters;
            labels(labels == jj) = idx;
            centroids(jj, :) = c_jj(1, :);
            centroids(num_clusters, :) = c_jj(2, :);
            if num_clusters == L
                break;
            end
            
        end
    end        
end

save('VQ.mat', 'centroids');
%{
figure(1);
for jj = 1:num_clusters
    plot(g_train(labels == jj, 1), g_train(labels == jj, 2), '.');
    hold on;
end
hold off;
%}
end