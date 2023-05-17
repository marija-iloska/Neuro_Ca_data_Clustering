clear all
close all
clc

% Load data
load Ca_data.mat

% Tastes in experiments
taste_strings = {'sugar', 'salt', 'inter'};

% Number of tastes in experiments
num_tastes = length(taste_strings);

% Number of trials PER taste
num_trials = [num_sugar_trials, num_salt_trials];

% Store frame indices for each trail
trial_frames = {sugar_start_end, salt_start_end};

% JOIN all "taste" trials
for i = 1 : num_tastes

    % Which taste
    trial = trial_frames{i};
    range = [];

    % Collect every trial index range
    for n = 1: num_trials(i)     
        temp = trial{n}(1) : trial{n}(2);
        range = [range temp];
    end

    % Get JOINT dFoF
    str = join(['joint_', taste_strings{i}]);
    assignin('base',str, dFoF(:, range));
end

