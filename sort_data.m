clear all
close all
clc

% Load data
load Ca_data.mat

% Tastes in experiments
intervals = {'sugar', 'salt'};

% Number of tastes in experiments
num_tastes = length(intervals);

% Number of trials PER taste
num_trials = [num_sugar_trials, num_salt_trials];

% Store frame indices for each trail
trial_frames = {sugar_start_end, salt_start_end};


% JOIN all "taste" trials
for i = 1 : num_tastes

    % Which taste
    trial = trial_frames{i};
    range = [];
    inter = [];

    % Collect every trial index range
    for n = 2: num_trials(i)     
        temp = trial{n}(1) : trial{n}(2);
        range = [range temp];
        temp_inter = trial{n-1}(end)+1 : trial{n}(1)-1; 
        inter = [inter, temp_inter];
        
    end

    % Get JOINT dFoF
    first_itrial = 1:trial{1}(1)-1;
    last_itrial = trial{num_trials(i)}(end)+1 :T;
    inter = [first_itrial, inter, last_itrial];  
    str = join(['inter_', intervals{i}]);
    assignin('base',str, inter);

    % Get JOINT dFoF
    range = [trial{1}(1) : trial{1}(2), range];
    str = join(['dFoF_', intervals{i}]);
    assignin('base',str, dFoF(:, range));
end

% add inter trial
intervals = {'sugar', 'salt', 'inter'};
frames_inter = intersect(inter_salt, inter_sugar);

str = join(['dFoF_', intervals{3}]);
assignin('base',str, dFoF(:, frames_inter));

save('trial_data.mat', 'dFoF_inter', 'dFoF_sugar', 'dFoF_salt', 'frames_inter', 'intervals')
