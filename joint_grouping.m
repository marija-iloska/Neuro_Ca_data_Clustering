clear all
close all
clc

% Load data
load Ca_data.mat
load trial_data.mat

% SUGAR, SALT, INTER
% Compute absolute area for each Neuron
L = round(0.25*length(dFoF_sugar(1,:)));
range1 = 2*L+1:3*L;
%range1 = 1500:3000;
for n = 1:num_neurons
    area_sugar(n) = trapz(abs(dFoF_sugar(n,range1)));
    area_salt(n) =  trapz(abs(dFoF_salt(n,range1)));
    area_inter(n) = trapz(abs(dFoF_inter(n,range1)));

end

% Sort area and get indices
[areaSU, idx_SU] = sort(area_sugar, 'descend');
[areaSA, idx_SA] = sort(area_salt, 'descend');
[areaIN, idx_IN] = sort(area_inter, 'descend');

% GET TWO CLUSTERS
dist_mat = pdist2(dFoF_salt, dFoF_salt);
Z = linkage(dist_mat, 'average');
figure
dendrogram(Z, 0)
set(gca, 'FontSize', 15)
idx = cluster(Z, 'Maxclust', 2);
top_salt = find(idx==1);

% Remove top neurons and recluster
FA_salt = dFoF_salt(setdiff(1:113, top_salt), :);

% GET TWO CLUSTERS
dist_mat = pdist2(FA_salt, FA_salt);
Z = linkage(dist_mat, 'average');
figure
dendrogram(Z, 0)
set(gca, 'FontSize', 15)
idx = cluster(Z, 'Maxclust', 2);
salt_new = find(idx==2);



% Get dF/F for top N most active
N = 20;
inters = intersect(idx_SU(1:N), idx_SA(1:N));
sep = inters;
%sep = setdiff(idx_SU(1:N), inters);

% dFoF of joint intervals
%F_SU = dFoF_sugar(idx_SU(sep), range1);
F_SA = dFoF_salt(idx_SA(sep), range1);
%F_IN = dFoF_inter(idx_IN(1:N), range1);


%% HIERACHICAL CLUSTERING
% Specify the number of clusters you want to obtain
k = 10;

% Perform hierarchical clustering using complete linkage
% link = { 'ward', 'complete', 'average', 'centroid', 'single'};
% L_SU = length(dFoF_sugar(1,:));
% L_IN = length(dFoF_inter(1,:));
% L_SA = length(dFoF_salt(1,:));
% range = 1:length(range1);
% for i = 1:1
%     Z_SA{i} = linkage(F_SA(:, range), link{i});
%     figure
%     dendrogram(Z_SA{i}, 0)
%     set(gca, 'FontSize', 15)
%     title(link{i}, 'FontSize', 20)
% end

% figure
% clust = [6, 8, 10];
% clust = top_sugar;
% for j = 1:length(clust)
%     plot(time(range), dFoF_sugar(clust(j),range), 'linewidth',1)
%     hold on
% end
% legend('3', '4',  'FontSize', 30)
% set(gca, 'FontSize', 20)

