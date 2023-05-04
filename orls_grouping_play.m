clear all
close all
clc

% Load data
load Ca_data.mat

range = 3000:6000;
idx = kmeans(dFoF(:, range), 10);

g = find(idx == 8);
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
xlim([550 750])
set(gca, 'FontSize', 30)
xlabel('Time', 'FontSize', 30)
ylabel('DeltaF/F', 'FontSize', 30)
