clear all
close all
clc

% Load data
load Ca_data.mat
%load trial_data_first.mat
%load trial_data_second.mat
load trial_data.mat

% FIRST for all
% SUGAR, SALT, INTER
%F = dFoF_first;
%F = dFoF_second;
F = {dFoF_sugar, dFoF_salt, dFoF_inter};
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
k = 4;

% GET TWO CLUSTERS
for m = 1 : M
    % Compute distance
    dist_mat = pdist2(F{m}(1:N,:), F{m}(1:N, :));
    %dist_mat = pdist(F{m}(1:N,:));

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

%     figure
%     dendrogram(Z, 0)
%     set(gca, 'FontSize', 15)
%     title(intervals{m}, 'FontSize', 15)
end


k = 4;
for m = 1:M
    dF = rel_dFoF{m}(:,:);

    % Compute distance
    dist_mat = pdist2(dF, dF);

    % Get linkage
    Z = linkage(dist_mat, 'ward');

    % figure
    % dendrogram(Z, 0)
    % set(gca, 'FontSize', 15)

    % Find cluster indices and store them
    idx_test = cluster(Z, 'maxclust', k);
    for kk = 1:k
        clusters{kk} = find(idx_test == kk);
    end

    % For each interval
    clust_int{m} = clusters;

    figure
    mm = 1;
    %mm = m;
    clust = datasample(clust_int{mm}, 3);
    dF = rel_dFoF{mm};
    clust = clust{1};
    range =1:1500;
    for j = 1:length(clust)
        plot(time(range), movmean(dF(clust(j),range),1), 'linewidth',2)
        hold on
    end
    legend('1', '2', 'FontSize', 30)
    title('Sugar, full interval', 'FontSize', 30)
    set(gca, 'FontSize', 20)


end





