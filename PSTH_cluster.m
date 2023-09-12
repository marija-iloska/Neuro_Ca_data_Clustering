clear all
close all
clc

% Load data
load Ca_data.mat
load neuron_PSTHs.mat
% load trial_data.mat
% 
% F_raw = {dFoF_sugar, dFoF_salt, dFoF_sugar, dFoF_salt};


intervals = {'Sucrose', 'Salt', 'Left', 'Right'};

% dFoF_salt = neurons(1:N, 1:50);
% dFoF_sugar = neurons(1:N, 51:100);
% dFoF_left = neurons(1:N, 101:150);
% dFoF_right= neurons(1:N, 151:200);
% F = {dFoF_sugar, dFoF_salt, dFoF_left, dFoF_right};

F = {neurons};

num_neurons = length(neurons(:,1));
N = num_neurons;
M = length(F);
% Compute absolute area for each Neuron
for m = 1:M

    corrF = corr(F{m}');

    for n = 1:num_neurons
        temp(n) = trapz(abs(F{m}(n,:)));
        corr_temp(n) = sum(abs(corrF(n,:)));
    end
    area{m} = temp;
    Fcor{m} = corr_temp;
    [area_sort{m}, idx_sort{m}] = sort(temp, 'descend');
    [Fcor_sort{m}, idx_sort_cor{m}] = sort(corr_temp, 'descend');

end


k = 2;

% GET TWO CLUSTERS
for m = 1 : M
    % Compute distance
    dist_mat = pdist2(F{m}(1:N,:), F{m}(1:N, :));


    % Get linkage
    Z = linkage(dist_mat, 'ward');

    % Find cluster indices
    idx_cluster{m} = cluster(Z, 'Maxclust', k);

    % Find cluster size and index of cluster size
    cluster_sz = [];
    for kk = 1:k
        cluster_sz = [cluster_sz, sum(idx_cluster{m} == kk)];
    end
    idx_big_group(m) = find(cluster_sz == max(cluster_sz));

    % Get the entire groups of neurons
    neurons_big_group{m} = find(idx_cluster{m} == idx_big_group(m));

    % Get intersections
    overlap{m} = intersect(neurons_big_group{m}, idx_sort_cor{m}(end -length(neurons_big_group{m}) : end));
    relevant{m} = setdiff(1:N, overlap{m});

    % Get rid of neurons
    %relevant{1} = 1:N;
    rel_dFoF{m} = F{m}(overlap{m}, :);

    figure
    dendrogram(Z, 0)
    set(gca, 'FontSize', 15)
    title(intervals{m}, 'FontSize', 15)
end

% SO FAR I JUST GOT RID OF SOME JUNK NEURONS




range = {1:50, 51:100, 101:150, 151:200, 1:100, 101:200, 1:200};

% Number of clusters
num_k = 2;
dF = rel_dFoF{1};

% Compute distance
dist_mat = pdist2(dF, dF);

% Get linkage
Z = linkage(dist_mat, 'ward');

% Find cluster indices and store them
idx_test = cluster(Z, 'maxclust', num_k);
for kk = 1:num_k
    clusters{kk} = find(idx_test == kk);
end


for n = 1:num_k

    % Cluster n 
    clust = clusters{n};

    indices = overlap{1}(clust)';
    subplot(2, 1, n)
    for j = 1:length(clust)
        plot(time(range{7}), movmean(dF(clust(j),range{7}),1), 'linewidth',2);
        hold on
    end
    plot(time(range{7}), mean(dF(clust, range{7}), 1), 'k', 'linewidth', 5);
    title(n, 'FontSize', 30)
    str = [string(indices), 'mean'];
    legend(str , 'FontSize', 15)
    set(gca, 'FontSize', 20)
end 

[~, idx_select] = sort([length(clusters{1}), length(clusters{2})]);
idx_select = idx_select(1);


r = 7;
dF_select = dF(clusters{idx_select}, range{r});

% Original indices
indices_select = overlap{1}(clusters{idx_select});


% Compute distance
dist_select = pdist2(dF_select, dF_select);

% Get linkage
Z_select = linkage(dist_select, 'ward');


idx_sel_clusters = cluster(Z_select, 'maxclust', 2);
for kk = 1:2
    clusters_select{kk} = find(idx_sel_clusters == kk);
end



for n = 1:2

    % Cluster n 
    clust = clusters_select{n};

    indices = indices_select(clust)';
    subplot(1,2, n)
    for j = 1:length(clust)
        plot(time(range{r}), movmean(dF_select(clust(j),range{r}),1), 'linewidth',2);
        hold on
    end
    plot(time(range{r}), mean(dF_select(clust, range{r}), 1), 'k', 'linewidth', 5);
    title(n, 'FontSize', 30)
    str = [string(indices), 'mean'];
    legend(str , 'FontSize', 15)
    set(gca, 'FontSize', 20)
end 






