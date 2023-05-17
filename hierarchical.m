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
[area, idx_sorted] = sort(areaN, 'descend');




% Get top N
N = num_neurons;

% N most active 
first = 1:N;

% N least active
last = num_neurons-N : num_neurons;


range = [];
for n = 1:10
    %temp = sugar_start_end{n}(1):sugar_start_end{n}(2);
    temp = salt_start_end{n}(1):salt_start_end{n}(2);
    range = [range temp];
end
range = 1:T;

F = dFoF(idx_sorted(first), range);

% Specify the number of clusters you want to obtain
k = 10;

% Compute the distance matrix between the rows of F
%dist_mat = pdist(F);
dist_mat = pdist2(F, F);


% Perform hierarchical clustering using complete linkage
Z = linkage(dist_mat, 'average');

% Extract the cluster assignments from the hierarchical tree
idx = cluster(Z, 'Maxclust', k);

% Plot a dendrogram of the hierarchical clustering result
figure(n)
dendrogram(Z, 0);
set(gca, 'FontSize', 15)

% % Display the cluster assignments
% disp(['Cluster assignments: ' num2str(idx')]);

figure(n+1)
Fplot = dFoF(idx_sorted, :);
clust = find(idx == 3);
%clust = [5, 8, 2, 12, 3, 7];
range = 2500:6500;
for j = 1:length(clust)
    plot(time(range), Fplot(clust(j),range), 'linewidth',1)
    hold on
end
pp = [149, 69, 247]/256;
legend('5', '8', '2', '12', 'FontSize', 20)

for n = 1:15
    xline(time(salt_start_end{n}(1)), 'r', 'LineWidth',3)
    hold on
    xline(time(salt_start_end{n}(2)), 'r', 'LineWidth',3)
    hold on
    xline(time(sugar_start_end{n}(1)), 'Color', pp, 'LineWidth',3)
    hold on
    xline(time(sugar_start_end{n}(2)), 'Color', pp, 'LineWidth',3)
    hold on
end
set(gca, 'FontSize', 20)
xlabel('Time', 'FontSize',20)
ylabel('dF/F', 'FontSize',20)
xlim([400, 900])



% NEW
group1 = find(idx == 2);
group2 = find(idx == 1);

% Active vs Inactive groups
F1 = Fplot(group1, :);
F2 = Fplot(group2, :);



dist_mat = pdist2(F2, F2);


% Perform hierarchical clustering using complete linkage
Z = linkage(dist_mat, 'average');

% Extract the cluster assignments from the hierarchical tree
idx = cluster(Z, 'Maxclust', k);

% Plot a dendrogram of the hierarchical clustering result
figure(n)
dendrogram(Z, 0);
set(gca, 'FontSize', 20)
title('Less Active', 'FontSize', 20)


figure(n+1)
%clust = find(idx == 3);
clust = [22, 4, 10];
range = 2500:6500;
for j = 1:length(clust)
    plot(time(range), F2(clust(j),range), 'linewidth',1)
    hold on
end
legend('22', '4', '10',  'FontSize', 30)
set(gca, 'FontSize', 20)


% for n = 1:15
%     xline(time(salt_start_end{n}(1)), 'r', 'LineWidth',3)
%     hold on
%     xline(time(salt_start_end{n}(2)), 'r', 'LineWidth',3)
%     hold on
%     xline(time(sugar_start_end{n}(1)), 'Color', pp, 'LineWidth',3)
%     hold on
%     xline(time(sugar_start_end{n}(2)), 'Color', pp, 'LineWidth',3)
%     hold on
% end
% pp = [149, 69, 247]/256;
% legend('8', '14', '12', 'FontSize', 20)
% xlim([400, 650])

