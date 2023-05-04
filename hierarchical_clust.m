clear all
close all
clc

load Ca_data.mat

% Specify the number of clusters you want to obtain
k = 10;


% Compute the distance matrix between the rows of F
range = 1:T;
F = dFoF(:,range);
dist_mat = pdist(F);
dist_mat = pdist2(F, F);
%dist_mat = 1- abs(corr(F'));

% Perform hierarchical clustering using complete linkage
Z = linkage(dist_mat, 'single');

% Extract the cluster assignments from the hierarchical tree
idx = cluster(Z, 'Maxclust', k);

% Plot a dendrogram of the hierarchical clustering result
dendrogram(Z);


i = datasample(1:k, 1);
g = find(idx == i);
range = 3000:5000;
more = range(end)+1 : range(end)+3000;
more =  [];
close all
% Plot all chosen neurons with high correlation to k1
for j = 1:length(g)
    plot(time([ range,more] ), dFoF(g(j),[range, more]), 'linewidth',1)
    hold on
end
pp = [149, 69, 247]/256;

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
xlim([450 650])
set(gca, 'FontSize', 30)
xlabel('Time', 'FontSize', 30)
ylabel('DeltaF/F', 'FontSize', 30)




figure(2)
range2 = 15000 : 17000;
for j = 1:length(g)
    plot(time([ range2,more] ), dFoF(g(j),[range2, more]), 'linewidth',1)
    hold on
end
pp = [149, 69, 247]/256;

for n = 30:45
    xline(time(salt_start_end{n}(1)), 'r', 'LineWidth',3)
    hold on
    xline(time(salt_start_end{n}(2)), 'r', 'LineWidth',3)
    hold on
    xline(time(sugar_start_end{n}(1)), 'Color', pp, 'LineWidth',3)
    hold on
    xline(time(sugar_start_end{n}(2)), 'Color', pp, 'LineWidth',3)
    hold on
end
xlim([1900 2200])
set(gca, 'FontSize', 30)
xlabel('Time', 'FontSize', 30)
ylabel('DeltaF/F', 'FontSize', 30)
