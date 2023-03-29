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
time = dt:dt:time_f;



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



% TOP DOWN BOTTOM UP BUT WITH GROUP AVERAGE ALSO
% Initialize
neuron_idx = 1:num_neurons;
area_Norg = area_N;
Acorr = corrcoef(area_N');
count = 1;
group_k = [];
while (isempty(neuron_idx) == 0)


    % Get idx of first top neuron
    if (length(neuron_idx) == num_neurons)
        k1 = 1;
    else
        k1 = datasample(neuron_idx, 1);
    end

    % Find max top 2 corr (self + one other)
    kmax = maxk(Acorr(k1,:), 2);
    
    % Get index of the other
    k_new = find( Acorr(k1,:) == kmax(end));


    if (ismember(k_new, group_k))
        ele = find(group_k == k_new);
        store_top{ele} = [store_top{ele}, k1];

        % Change in area
        area_N(ele,:) = mean(area_N(store_top{ele}, :), 1);
    else
        % Store Group
        store_top{count} = [k1, k_new];
        group_k = [group_k, k1];

        % Change in area
        area_N(k1,:) = mean(area_N(store_top{count}, :), 1);
        % Increase count
        count = count + 1;
    end

    % Remove area for selected element
    area_N(k_new,:) = NaN;

    % Remove index from exploration space
    neuron_idx = setdiff(neuron_idx, [k1, k_new]);

    % Recompute Corr
    Acorr = corrcoef(area_N');

end

% Visualize
% sample random group
r = datasample(1:10, 1);
%idx = r*1000:(r*1000 + 3000);
idx = 1:2000;
j = datasample(1:length(store_top), 1);
for k = 1:length(store_top{j})
    figure(1)
    plot(area_Norg(store_top{j}(k), 2000:3000))
    hold on
%     figure(2)
%     plot(dFoF(store_top{j}(k), idx))
%     hold on
end
plot(mean(area_Norg(store_top{j}, 2000:3000), 1), 'k', 'linewidth',2)





