clear all
close all
clc

% Load data
load Ca_data.mat
load trial_data.mat

% SUGAR, SALT, INTER
% Compute absolute area for each Neuron
for n = 1:num_neurons
    area_sugar(n) = trapz(abs(dFoF_sugar(n,:)));
    area_salt(n) =  trapz(abs(dFoF_salt(n,:)));
    area_inter(n) =  trapz(abs(dFoF_inter(n,:)));

end

% Sort area and get indices
[areaSU, idx_SU] = sort(area_sugar, 'descend');
[areaSA, idx_SA] = sort(area_salt, 'descend');
[areaIN, idx_IN] = sort(area_inter, 'descend');

% Get dF/F for top N most active
N = 10;
F_SU = dFoF_sugar(idx_SU(1:N), :);
F_SA = dFoF_salt(idx_SA(1:N), :);
F_IN = dFoF_inter(idx_IN(1:N), :);


% Specify the number of clusters you want to obtain
k = 10;

% Compute the distance matrix between the rows of F
%dist_mat = pdist(F);
dist_mat = pdist2(F_SU, F_SU);

% Perform hierarchical clustering using complete linkage
Z = linkage(dist_mat, 'average');

% Extract the cluster assignments from the hierarchical tree
idx = cluster(Z, 'Maxclust', k);

% Plot a dendrogram of the hierarchical clustering result
figure(n)
dendrogram(Z, 0);
set(gca, 'FontSize', 15)

