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
T = length(dFoF(1,:));
time_f = T*dt;

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
j = datasample(1:113, 8);
idx = 1:3000;
for k = 1:8
    subplot(4, 2, k)
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



% K-means clustering based on area
num_clusters = 20;
ktest = kmeans([N_area', ones(1, num_neurons)'], num_clusters);


% ABSOLUTE AREA PLOT
for k = 1:num_clusters
    id{k} = find(ktest==k);
end
plot(N_area, '.', 'Linewidth', 2)
hold on
for k = 1:num_clusters
    yline(N_area(id{k}(1)))
    hold on
end
set(gca,'FontSize', fsz)
xlabel('Neuron','FontSize', fsz)
ylabel('Absolute Area', 'FontSize', fsz)
legend('Sorted order','FontSize', fsz)


idx = 20000:21000;
figure(5)
for k = 1:num_clusters
    plot(time(idx), mean(dFoF(id{k},idx), 1))
    hold on
end




% Plot bottom 4 based on area
S = 10;
idx = 20000:21500;
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
idx = 20000:21500;
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


% Plot bottom 4 based on area
S = 10;
idx = 20000:21500;
idx_new = [idx_sort(1:3), idx_sort(end-3:end)];
figure(4)
for i = 1 : 6
    subplot(3, 2, i)
    plot(time(idx), dFoF(idx_new(i),idx), 'k')
    set(gca,'FontSize', fsz)
    xlabel('Time','FontSize', fsz)
    ylabel('DeltaF/F', 'FontSize', fsz)
    ylim([-0.5, 4])
    xlim([time(idx(1)), time(idx(end))])
end


% ft = fft(dFoF);
% 
% corrcoef(ft');
% % FFT test
% ft{1} = fft(dFoF(idx_sort(1), :));
% ft{2} = fft(dFoF(50, :));
% ft{3} = fft(dFoF(idx_sort(3), :));
% ft{4} = fft(dFoF(idx_sort(end-1), :));
% j = [1,113];
% for k = 1:2
%     plot(time(idx), movmean(dFoF(idx_sort(j(k)),idx),1))
%     hold on
% end
% 



% Get correlation
store = [];
nn = 4;
n = 50;
test = corrcoef(dFoF(idx_sort,:)');
kmax = maxk(test(n,:), nn);
for k = 1:nn
    store(k) = find(test(n,:) == kmax(k));
end

idx = 4000:6000;
for k = 1:nn
    subplot(round(nn/2),2,k)
    plot(dFoF(store(k), idx))
    hold on
end


idx = 1:3500;
S = 10;
figure(3)
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




% ABSOLUTE AREA in TIME STEPS
num_groups = 3;
n=1;
i = 0.5;
step = i;
count = 1;

clear area_N
% For each time step
for t = 1:T

    % For each neuron, get number of spikes within time bin

    % Find max rate
    if (time(t) > i)
        
        for j = 1:num_neurons
        
            area_N(j,count) = trapz(abs(dFoF(j,n:t)));
        end
        i = i + step;
        n = t;
        count = count+1;
    end
    
end



test = corrcoef(area_N');
% Get correlation
store = [];
nn = 2;
n = 113;
kmax = maxk(test(n,:), nn);
for k = 1:nn
    store(k) = find(test(n,:) == kmax(k));
end

idx = 17000:20000;
for k = 1:nn
    figure(6)
    plot(area_N(store(k),1000:3000))
    hold on

    figure(7)
    plot(time(idx), dFoF(store(k),idx))
    hold on
end

% TOP DOWN BOTTOM UP
% Initialize
neuron_idx = 1:num_neurons;
Acorr = corrcoef(area_N');
count = 1;
while (isempty(neuron_idx) == 0)

    % Get index for top neuron
    k1 = neuron_idx(1);

    % Max corrs for top and bottom neurons based on area
    kmax_top = maxk(Acorr(k1,:), 2);
    store_top{count} = [k1, find(Acorr(k1,:) == kmax_top(end))];
    neuron_idx = setdiff(neuron_idx, store_top{count});

    % Get index for bottom neurons
    k2 = neuron_idx(end);

    kmax_bot = maxk(Acorr(k2,:), 2);
    store_bot{count} = [k2, find(Acorr(k2,:) == kmax_bot(end))];
    neuron_idx = setdiff(neuron_idx, store_bot{count});
    
    count = count + 1;

    if( length(neuron_idx) == 1)
        break
    end

end

% TEST
j = datasample(1:42, 1);
idx = 17000:20000;
for k = 1:2
    figure(6)
    plot(area_N(store_bot{j}(k), 1000:3000))
    hold on

    figure(7)
    plot(time(idx), dFoF(store_bot{j}(k),idx))
    hold on
end


