clear all
close all
clc

% Load data
load Ca_data.mat


%----------------------------------%
%    GET ABSOLUTE AREA in steps
%----------------------------------%

% Initialize
n=1;    count = 1;
i = 0.3;    step = i;

clear area_N
while (i <= time_f)

    % Get time index
    t = find(time > i);
    t = t(1);

    % Integrate - Absolute Area for every Neuron with steps    
    area_N(:, count) = trapz(abs(dFoF(:,n:t)'));

    % Update indices
    i = i + step;
    n = t+1;
    count = count + 1;
 
end





%----------------------------------%
%        AREA grouping
%----------------------------------%
[store] = grouping(area_N, num_neurons);

% PLOTTING grouping based on AREA
pp = [149, 69, 247]/256;
idx = 1:2500;

% Sample random pair
j = datasample(1:length(store), 1);
for k = 2:length(store{j})
%     figure(1)
%     plot(area_Norg(store_top{j}(k), :))
%     hold on
    figure(2)
    plot(time(idx), movmean(dFoF(store{j}(k), idx),1))
    hold on
end
% for n = 1:4
%     figure(2)
%     xline(time(salt_start_end{n}(1)), 'k', 'LineWidth',2)
%     hold on
%     xline(time(salt_start_end{n}(2)), 'k', 'LineWidth',2)
%     hold on
%     xline(time(sugar_start_end{n}(1)), 'Color', pp, 'LineWidth',2)
%     hold on
%     xline(time(sugar_start_end{n}(2)), 'Color', pp, 'LineWidth',2)
%     hold on
% end
xlabel('Time', 'FontSize', 20)
ylabel('DeltaF/F', 'FontSize', 20)
ylabel('Absolute Area', 'FontSize', 20)





%----------------------------------%
%        1 NEURON grouping
%----------------------------------%

% Choose length of interval and step size
i = 80;     step = i;

% Choose neuron of interest
k1 = 16;
str = join(['Neuron ', num2str(k1)]);

% Initialize
count = 1;  n = 1;  group = [];
clear time_idx


% Group k1 with highest correlation within intervals
while (i < time_f)

    % Get time index
    t = find(time > i);
    t = t(1);

    % Integrate - Absolute Area for every Neuron with steps    
    Fcorr = corrcoef(dFoF(:,n:t)');

    kmax = maxk(Fcorr(k1,:), 2);
    k_new = find(Fcorr(k1,:) == kmax(end));
    group{count} = [k1, k_new];

    % Update indices
    count = count + 1;
    i = i + step;
    time_idx(count) = t;
    n = t+1;
end


% PLOTTING FOR 1 NEURON GROUPING
n = 1;  count = 1;

% Moving average degree
deg = 3;

% Plot all chosen neurons with high correlation to k1
for j = 1:15
    figure(5)
    n_idx = group{j}(2);
    t = time_idx(count+1);
    plot(time(n:t), movmean(dFoF(n_idx,n:t), deg), 'linewidth',2)
    hold on
    n = time_idx(count+1)+1;
    count = count + 1;
end

% Plot k1
figure(5);
plot(time(1:t), movmean(dFoF(k1,1:t), deg), 'k', 'linewidth',2)
hold on

% Plot Trial intervals
for n = 1:15
    figure(5)
    xline(time(salt_start_end{n}(1)), 'r', 'LineWidth',2)
    hold on
    xline(time(salt_start_end{n}(2)), 'r', 'LineWidth',2)
    hold on
    xline(time(sugar_start_end{n}(1)), 'Color', pp, 'LineWidth',2)
    hold on
    xline(time(sugar_start_end{n}(2)), 'Color', pp, 'LineWidth',2)
    hold on
end

% Aesthetics
%xlim([500, 800])
xlim([0 500])
set(gca, 'FontSize', 30)
xlabel('Time', 'FontSize', 30)
ylabel('DeltaF/F', 'FontSize', 30)
title(str, 'FontSize', 30)


