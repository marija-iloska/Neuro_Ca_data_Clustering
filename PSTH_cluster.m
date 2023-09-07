clear all
close all
clc

% Load data
load Ca_data.mat
load neuron_PSTHs.mat
% load trial_data.mat
% 
% F_raw = {dFoF_sugar, dFoF_salt, dFoF_sugar, dFoF_salt};

N = num_neurons;

intervals = {'Sucrose', 'Salt', 'Left', 'Right'};

dFoF_salt = neurons(1:N, 1:50);
dFoF_sugar = neurons(1:N, 51:100);
dFoF_left = neurons(1:N, 101:150);
dFoF_right= neurons(1:N, 151:200);
F = {dFoF_sugar, dFoF_salt, dFoF_left, dFoF_right};


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
    rel_dFoF{m} = F{m}(relevant{m}, :);

%     figure
%     dendrogram(Z, 0)
%     set(gca, 'FontSize', 15)
%     title(intervals{m}, 'FontSize', 15)
end





num_k = [5,5,5,5];
for m = 1:M
    dF = rel_dFoF{m};

    % Compute distance
    dist_mat = pdist2(dF, dF);

    % Get linkage
    Z = linkage(dist_mat, 'ward');

    % Find cluster indices and store them
    idx_test = cluster(Z, 'maxclust', num_k(m));
    for kk = 1:num_k(m)
        clusters{kk} = find(idx_test == kk);
    end

    % For each interval
    clust_int{m} = clusters;

    figure
    mm = m;
    clust = datasample(clust_int{mm}, 3);
    dF = rel_dFoF{mm};
    clust = clust{1};
    % Tracking indices
    indices = relevant{m}(clust);
    range = 1:50;
    for j = 1:length(clust)
        plot(time(range), movmean(dF(clust(j),range),1), 'linewidth',2);
        hold on
    end
    plot(time(range), mean(dF(clust, range), 1), 'k', 'linewidth', 5);
    title(intervals{mm}, 'FontSize', 30)
    str = [string(indices), 'mean'];
    legend(str , 'FontSize', 15)
    set(gca, 'FontSize', 20)

    filename = join(['figs/', intervals{mm}, '.eps']);
    print(gcf, filename, '-depsc2', '-r300');

end





