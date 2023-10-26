function [Z_store, overlap, rest] = clustering(F, N, M, k, link, idx_sort)

% GET TWO CLUSTERS
for m = 1 : M
    % Compute distance
    dist_mat = pdist2(F, F(1:N, :));
    %dist_mat = pdist(F{m}(1:num_neurons,:));

    % Get linkage
    Z = linkage(dist_mat, link);
    Z_store{m} = Z;

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
    rest{m} = setdiff(1:N, overlap{m});


end




end