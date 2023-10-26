clear all
close all
clc

% Choose session
session = {'s202', 's301', 's302', 's313', 'sf203', 'sf309', 'sf311'};

idx = 2;
ca_str = join(['Ca_Data/Ca_', session{idx}, '.mat']);
tr_str = join(['Sort_Data/trial_', session{idx}, '.mat']);

% Load data
load(ca_str)
load(tr_str)

% For my PSTHs a good clustering linkage is 'ward'
% PSTHs concatenate as [sugar, salt, left, right]

F = {psth_taste{1}, psth_taste{2}, psth_decis{1}, psth_decis{2}};
neurons = [psth_taste{1}, psth_taste{2}, psth_decis{1}, psth_decis{2}];
%load neuron_PSTHs.mat

% Relevant Sections
intervals = {'Sucrose', 'Salt', 'Left', 'Right'};
M = length(intervals);


%% ABSOULTE AREA

% Compute absolute area for each Neuron
for m = 1:M
    for n = 1:num_neurons
        temp(n) = trapz(abs(F{m}(n,:)));
    end
    area{m} = temp;
    [area_sort{m}, idx_sort{m}] = sort(temp, 'descend');
end



%% GET RID OF LOUD MAJORITY
% Linkage
link = 'average';

% Number of clusters
k = 3;

% Overall clustering and overlap with area
[Z, overlap, rest] = clustering(F, num_neurons, M, k, link, idx_sort);

% Visualize
figure
m = 3;
dendrogram(Z{m}, 0)
set(gca, 'FontSize', 15)
title(intervals{m}, 'FontSize', 15)



%% CLUSTERING ON SMALL GROUP

% Original indices would be overlap{m}(clusters{k1, k2, ..})

% Number of clusters to break in
num_k = 20;

% Choose index [sugar, salt, left, right]
m = 4;

% Only Interval data here
%dF = neurons(:, 136:180);
dF = F{m}(overlap{m}, :);

% Compute distance
dist_mat = pdist2(dF, dF);

% Get linkage
Z = linkage(dist_mat, 'average');

% figure
% dendrogram(Z, 0)
% set(gca, 'FontSize', 15)

% Find cluster indices and store them
idx_test = cluster(Z, 'maxclust', num_k);
for kk = 1:num_k
    clusters{kk} = find(idx_test == kk);
end


%% SAMPLE A CLUSTER AND PLOT

figure;
for kk = 1:num_k
    % Sample one of the clusters
    %clust = datasample(clusters, 1);
    % Extract the indices in the cluster
    %clust = clust{1};

    clust = clusters{kk};

    % Tracking indices
    %indices = overlap{m}(clust)';
    indices = clust';


    subplot(4, 5, kk)
    for j = 1:length(clust)
        plot(neurons(indices(j),:), 'linewidth',2);
        hold on
    end
    plot(mean(neurons(indices, :), 1), 'k', 'linewidth', 2);
    %title(intervals{m}, 'FontSize', 30)
    str = [string(indices), 'mean'];
    legend(str , 'FontSize', 15)
    set(gca, 'FontSize', 20)


end
% filename = join(['figs/', intervals{mm}, '.eps']);
% print(gcf, filename, '-depsc2', '-r300');






