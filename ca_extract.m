clear all
close all
clc

% Choose which data to load
choice = {'neurons_aligned_FL', 'neurons_aligned_FC', 'neurons_aligned_valve'};
idx = 2;

% Load and convert data
neuron_data = struct2cell(load(choice{idx}));
neuron_data = neuron_data{:};


% Number of Neurons
num_neurons = length(neuron_data);

% Read in dFoF for each neuron
for n = 1:num_neurons
   dFoF(n,:) = neuron_data(n).F;
end

% Time steps
dt = 0.129;

% Length of time series
T = length(dFoF(1,:));
time_f = T*dt;

% Create time array
time = dt:dt:T;


% Number of trials where Mouse chose Salt vs Sugar
num_salt_trials = length(neuron_data(1).NaCl.spikes);
num_sugar_trials = length(neuron_data(1).Sucrose.spikes);


% Get trial start and ends for NaCl 
for n = 1:num_salt_trials
    temp = neuron_data(1).NaCl.Frames_index{n};
    salt_start_end{n} = [temp(1), temp(end)];
end

% Get trial start and ends for Sucrose
for n = 1:num_sugar_trials
    temp = neuron_data(1).Sucrose.Frames_index{n};
    sugar_start_end{n} = [temp(1), temp(end)];
end


% Save CA data
save("Ca_data.mat", 'dt', 'T', 'time', 'time_f', 'dFoF', 'salt_start_end', 'salt_start_end', 'num_neurons')






% pp = [149, 69, 247]/256;
% tq = [65, 150, 186]/256;

% Visualization
% fsz = 15;
% j = datasample(1:113, 8);
% idx = 1:3000;
% for k = 1:8
%     subplot(4, 2, k)
%     plot(time(idx), dFoF(j(k),idx), 'k')
%     hold on
%     for n = 1:3
%         xline(time(salt_start_end{n}(1)), 'r', 'LineWidth',2)
%         hold on
%         xline(time(salt_start_end{n}(2)), 'r', 'LineWidth',2)
%         hold on
%         xline(time(sugar_start_end{n}(1)), 'Color', pp, 'LineWidth',2)
%         hold on
%         xline(time(sugar_start_end{n}(2)), 'Color', pp, 'LineWidth',2)
%         hold on
%     end
%     set(gca, 'FontSize', fsz)
%     xlabel('Time','FontSize', fsz)
%     ylabel('DeltaF/F','FontSize', fsz)
% end












