clear all
close all
clc

% Load data
load Ca_data.mat
%load trial_data_first.mat
load trial_data_second.mat

% FIRST for all
% SUGAR, SALT, INTER 
F = dFoF_second;
M = length(F);
% Compute absolute area for each Neuron
for m = 1:M
    for n = 1:num_neurons
        temp(n) = trapz(abs(F{m}(n,:))); 
    end
    area{m} = temp;
    [area_sort{m}, idx_sort{m}] = sort(temp, 'descend');
end
N = 113;

% GET TWO CLUSTERS
for m = 1 : M
    % Compute distance
    dist_mat = pdist2(F{m}(1:N,:), F{m}(1:N, :));

    % Get linkage
    Z = linkage(dist_mat, 'average');
    
    % Find cluster indices
    idx_cluster{m} = cluster(Z, 'Maxclust', 2);

    % Find cluster size and index of cluster size
    cluster_sz = [sum(idx_cluster{m} == 1), sum(idx_cluster{m}==2)];
    idx_small_group(m) = find(cluster_sz == min(cluster_sz));

    % Get the entire groups of neurons
    neurons_small_group{m} = find(idx_cluster{m} == idx_small_group(m));

    % Get intersections
    overlap{m} = intersect(neurons_small_group{m}, idx_sort{m}(1:length(neurons_small_group{m})));

    figure
    dendrogram(Z, 0)
    set(gca, 'FontSize', 15)
    title(intervals{m}, 'FontSize', 15)
end




% VISUALIZE AREA
for m = 1:length(F)
    figure(1)
    plot(area{m}, 'linewidth',1)
    hold on

    figure(2)
    plot(area_sort{m}, 'linewidth',1)
    hold on


end
xline(length(neurons_small_group{1}), 'color', 'k', 'linewidth',1)
hold on
xline(length(neurons_small_group{2}), 'color', 'b', 'linewidth',1)
hold on
xline(length(neurons_small_group{3}), 'color', 'm', 'linewidth',1)
hold on
legend('sugar', 'salt', 'inter', 'FontSize', 15)





