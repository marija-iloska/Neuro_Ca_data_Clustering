clear all
close all
clc

% Load data
load Ca_data.mat

% Tastes in experiments
intervals = {'sugar', 'salt'};
split = {'entire', 'first', 'second'};


% Number of tastes in experiments
num_tastes = length(intervals);

% Number of trials PER taste
num_trials = [num_sugar_trials, num_salt_trials];

% Store frame indices for each trail
trial_frames = {sugar_start_end, salt_start_end};

% Normalizing sum
min_F = min(min(dFoF));
w_norm = sum(dFoF - min_F, 2);
dFoF_norm = (dFoF -min_F)./w_norm;
dFoF_norm = dFoF_norm/(max(max(dFoF_norm)));


% JOIN all "taste" trials
for i = 1 : num_tastes

    % Which taste
    trial = trial_frames{i};
    range = [];
    first = [];
    second = [];
    psth_entire = [];
    psth_first  = [];
    psth_second = [];

    % Collect every trial index range - ENTIRE trial range
    for n = 1: num_trials(i)   

        % Extract range
        temp =  trial{n}(1): trial{n}(2);
        first_temp = temp(1) : floor(0.5*(temp(1) + temp(end)));
        second_temp = setdiff(temp, first_temp);
        

        % Concatenate with all previous extracted ranges
        range = [range temp];  
        first = [first first_temp];
        second = [second second_temp];
        
        % PSTH
        psth_entire = [psth_entire;   dFoF_norm(:,temp)];
        psth_first  = [psth_first;    dFoF_norm(:,first_temp)];
        psth_second = [psth_second;   dFoF_norm(:,second_temp)];
    end

    % Get JOINT dFoF ranges
    entire_trial{i} = range;
    first_half{i} = first;
    second_half{i} = second;

    % PSTH
    psth_trial{i} = squeeze(mean(psth_entire, 1));
    psth_taste{i} = squeeze(mean(psth_first, 1));
    psth_decis{i} = squeeze(mean(psth_second, 1));


    % Name and get values using indices
    str = join(['dFoF_entire', '_', intervals{i}]);
    assignin('base',str, dFoF(:, entire_trial{i}));

    str = join(['dFoF_first', '_', intervals{i}]);
    assignin('base',str, dFoF(:, first_half{i}));

    str = join(['dFoF_second', '_', intervals{i}]);
    assignin('base',str, dFoF(:, second_half{i}));

end


% Get joint inter trial frames
inter_trial = setdiff(1:T, [entire_trial{1}, entire_trial{2}]);

str = 'dFoF_inter';
assignin('base',str, dFoF(:, inter_trial));


% Save variables
frames = {entire_trial, first_half, second_half, inter_trial};
dFoF_sugar = {dFoF_entire_sugar, dFoF_first_sugar, dFoF_second_sugar};
dFoF_salt = {dFoF_entire_salt, dFoF_first_salt, dFoF_second_salt};

% Save file
save('trial_data.mat', 'dFoF_sugar', 'dFoF_salt', 'frames', 'split', ...
    'dFoF_norm', 'psth_decis', 'psth_taste', 'psth_trial')


