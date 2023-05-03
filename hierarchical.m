clear all
close all
clc

% Load data
load Ca_data.mat

% Compute absolute area for each Neuron
for n = 1:num_neurons
    areaN(n) = trapz(abs(dFoF(n,:)));
end

% Sort area and get indices
[area, idx] = sort(areaN, 'descend');

% Get top 50
F = dFoF(idx(1:40), :);

% Specify the number of clusters you want to obtain
k = 3;

% Compute the distance matrix between the rows of F
dist_mat = pdist(F);

% Perform hierarchical clustering using complete linkage
Z = linkage(dist_mat, 'complete');

% Extract the cluster assignments from the hierarchical tree
idx = cluster(Z, 'Maxclust', k);

% Plot a dendrogram of the hierarchical clustering result
dendrogram(Z, 0);

% Display the cluster assignments
disp(['Cluster assignments: ' num2str(idx')]);