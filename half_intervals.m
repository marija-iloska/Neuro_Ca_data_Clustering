clear all
close all
clc

% Load data
load Ca_data.mat
load trial_data_first.mat
%load trial_data_second.mat

% FIRST for all
% SUGAR, SALT, INTER 
%F = dFoF_first;
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
k = 10;

% GET TWO CLUSTERS
for m = 1 : M
    % Compute distance
    dist_mat = pdist2(F{m}(1:N,:), F{m}(1:N, :));

    % Get linkage
    Z = linkage(dist_mat, 'average');
    
    % Find cluster indices
    idx_cluster{m} = cluster(Z, 'Maxclust', k);

    % Find cluster size and index of cluster size
    cluster_sz = [];
    for kk = 1:k
        cluster_sz = [cluster_sz, sum(idx_cluster{m} == kk)];
    end
    idx_big_group(m) = find(cluster_sz == max(cluster_sz));
    %idx_small_group(m) = find(cluster_sz == min( setdiff(cluster_sz, min(cluster_sz))));

    % Get the entire groups of neurons
    neurons_big_group{m} = find(idx_cluster{m} == idx_big_group(m));
    %neurons_small_group{m} = find(idx_cluster{m} == idx_small_group(m));

    % Get intersections
    %overlap{m} = intersect(neurons_small_group{m}, idx_sort{m}(1:length(neurons_small_group{m})));
    overlap{m} = intersect(neurons_big_group{m}, idx_sort{m}(end -length(neurons_big_group{m}) : end));
    relevant{m} = setdiff(1:N, overlap{m});

    % Get rid of neurons
    rel_dFoF{m} = F{m}(relevant{m}, :);

    figure
    dendrogram(Z, 0)
    set(gca, 'FontSize', 15)
    title(intervals{m}, 'FontSize', 15)
end




F = rel_dFoF{1}(:,1:600);

% Compute distance
dist_mat = pdist2(F, F);

% Get linkage
Z = linkage(dist_mat, 'ward');

figure
dendrogram(Z, 0)
set(gca, 'FontSize', 15)


figure(n+1)
%clust = find(idx == 3);
clust = [13, 41, 5, 10, 25, 33 ];
%clust = [4, 21, 35, 36, 37];
clust = [11, 14, 15, 45];
% clust = [12, 32, 31, 19];
% clust = [24, 44, 26, 38]; %, 2, 26, 38];
clust = [6, 39, 27];
range = 1 : 500;
for j = 1:length(clust)
    plot(time(range), movmean(F(clust(j),range),1), 'linewidth',2)
    hold on
end
%legend('22', '4', '10',  'FontSize', 30)
set(gca, 'FontSize', 20)



% VISUALIZE AREA
% for m = 1:length(F)
%     figure(1)
%     plot(area{m}, 'linewidth',1)
%     hold on
% 
%     figure(2)
%     plot(area_sort{m}, 'linewidth',1)
%     hold on
% 
% 
% end
% xline(length(neurons_big_group{1}), 'color', 'k', 'linewidth',1)
% hold on
% xline(length(neurons_big_group{2}), 'color', 'b', 'linewidth',1)
% hold on
% xline(length(neurons_big_group{3}), 'color', 'm', 'linewidth',1)
% hold on
% legend('sugar', 'salt', 'inter', 'FontSize', 15)





