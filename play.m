clear all
close all
clc

% Load data
load neurons_aligned_FC.mat

% Read in F
num_neurons = length(neurons_aligned_FC);

for n = 1:num_neurons
   dFoF(n,:) = neurons_aligned_FC(n).F;
end

% Time steps
dt = 0.129;

% Length of time series
T = length(dFoF(1,:))*dt;

% Create time array
time = dt:dt:T;


% Number of trials where Mouse chose Salt vs Sugar
num_salt_trials = length(neurons_aligned_FC(1).NaCl.spikes);
num_sugar_trials = length(neurons_aligned_FC(1).Sucrose.spikes);

% Get trial start and ends for NaCl 
for n = 1:num_salt_trials
    temp = neurons_aligned_FC(1).NaCl.Frames_index{n};
    salt_start_end{n} = [temp(1), temp(end)];
end

% Get trial start and ends for Sucrose
for n = 1:num_sugar_trials
    temp = neurons_aligned_FC(1).Sucrose.Frames_index{n};
    sugar_start_end{n} = [temp(1), temp(end)];
end


pp = [149, 69, 247]/256;
tq = [65, 150, 186]/256;

% Visualization
fsz = 15;
j = datasample(1:113, 4);
idx = 1:3000;
for k = 1:4
    subplot(2, 2, k)
    plot(time(idx), dFoF(j(k),idx), 'k')
    hold on
    for n = 1:3
        xline(time(salt_start_end{n}(1)), 'r', 'LineWidth',2)
        hold on
        xline(time(salt_start_end{n}(2)), 'r', 'LineWidth',2)
        hold on
        xline(time(sugar_start_end{n}(1)), 'Color', pp, 'LineWidth',2)
        hold on
        xline(time(sugar_start_end{n}(2)), 'Color', pp, 'LineWidth',2)
        hold on
    end
    set(gca, 'FontSize', fsz)
    xlabel('Time','FontSize', fsz)
    ylabel('DeltaF/F','FontSize', fsz)
end




% Find max and min and mean and absolute area in every neuron
for n = 1:num_neurons
    
    n_max(n) = max(dFoF(n,:));
    n_min(n) = min(dFoF(n,:));
    n_mean(n) = mean(dFoF(n,:));
    n_std(n) = std(dFoF(n,:));
    n_area(n) = trapz(abs(dFoF(n,:)));

end

% Sort area highest to lowest
[N_area, idx_sort] = sort(n_area, 'descend');


% Basic statistics PLOT
plot(n_mean, 'k', 'LineWidth',2)
hold on
plot(n_max, 'Linewidth', 2)
hold on
plot(n_min, 'LineWidth',2)
set(gca,'FontSize', fsz)
xlabel('Neuron','FontSize', fsz)
ylabel('Intensity', 'FontSize', fsz)
legend('Mean', 'Max', 'Min','FontSize', fsz)


% ABSOLUTE AREA PLOT
plot(N_area, 'Linewidth', 2)
set(gca,'FontSize', fsz)
xlabel('Neuron','FontSize', fsz)
ylabel('Absolute Area', 'FontSize', fsz)
legend('Sorted order','FontSize', fsz)



% Plot bottom 4 based on area
S = 10;
idx = 7500:9500;
figure(1)
for i = 1 :S
    subplot(round(S/2), 2, i)
    plot(time(idx), dFoF(idx_sort(end-i-S),idx), 'k')
    set(gca,'FontSize', fsz)
    xlabel('Time','FontSize', fsz)
    ylabel('DeltaF/F', 'FontSize', fsz)
    ylim([-0.5, 4])
    xlim([time(idx(1)), time(idx(end))])
end


% Plot bottom 4 based on area
S = 10;
idx = 7500:9500;
figure(2)
for i = 1:S
    subplot(round(S/2), 2, i)
    plot(time(idx), dFoF(idx_sort(i),idx), 'Color', tq)
    set(gca,'FontSize', fsz)
    xlabel('Time','FontSize', fsz)
    ylabel('DeltaF/F', 'FontSize', fsz)
    ylim([-0.5, 4])
    xlim([time(idx(1)), time(idx(end))])
end



idx = 1:3500;
S = 30;
plot(time(idx), mean(dFoF(idx_sort(1:S),idx), 1))
hold on
plot(time(idx), mean(dFoF(idx_sort(end-S:end),idx), 1), 'k')
hold on
for n = 1:7
    xline(time(salt_start_end{n}(1)), 'r', 'LineWidth',2)
    hold on
    xline(time(salt_start_end{n}(2)), 'r', 'LineWidth',2)
    hold on
%     xline(time(sugar_start_end{n}(1)), 'Color', pp, 'LineWidth',2)
%     hold on
%     xline(time(sugar_start_end{n}(2)), 'Color', pp, 'LineWidth',2)
%     hold on
end

