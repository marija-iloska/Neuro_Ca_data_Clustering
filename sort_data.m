clear all
close all
clc

% Load data
load Ca_data.mat

% Tastes in experiments
intervals = {'sugar', 'salt'};
split = {'first', 'second'};
choice = split{1};

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
        temp =  (trial{n}(1) + floor( 0.5*(trial{n}(2) - trial{n}(1)))) : trial{n}(2);
        %temp = trial{n}(1) : (trial{n}(1) + floor( 0.5*(trial{n}(2) - trial{n}(1))))
        %temp = trial{n}(1) : trial{n}(2);
        range = [range temp];
        % temp_inter = trial{n-1}(end)+1 : trial{n}(1)-1;
        temp_inter = (trial{n-1}(end)+1 + floor( 0.5*(trial{n}(1)-1 - trial{n-1}(end)+1) )) : trial{n}(1)-1; 
        inter = [inter, temp_inter];
        
    end

    % Get JOINT dFoF
%     first_itrial = 1:trial{1}(1)-1;
%     last_itrial = trial{num_trials(i)}(end)+1 :T;
%     inter = [first_itrial, inter, last_itrial];  
    str = join([ choice, '_inter_', intervals{i}]);
    assignin('base',str, inter);

    % Get JOINT dFoF
    %range = [trial{1}(1) : trial{1}(2), range];
    range = [trial{1}(1) : (trial{1}(1) + floor(0.5*(trial{1}(2) - trial{1}(1)))), range];
    str = join(['dFoF_', choice, '_', intervals{i}]);
    assignin('base',str, dFoF(:, range));
end

% add inter trial
intervals = {'sugar', 'salt', 'inter'};
frames_inter = intersect(first_inter_salt, first_inter_sugar);
%frames_inter = intersect(second_inter_salt, second_inter_sugar);

str = join(['dFoF_', choice, '_', intervals{3}]);
assignin('base',str, dFoF(:, frames_inter));

dFoF_first = {dFoF_first_sugar, dFoF_first_salt, dFoF_first_inter};
%dFoF_second = {dFoF_second_sugar,dFoF_second_salt, dFoF_second_inter};


%save('trial_data.mat', 'dFoF_inter', 'dFoF_sugar', 'dFoF_salt', 'frames_inter', 'intervals')
%save('trial_data_first.mat', 'dFoF_first', 'dFoF_first_inter', 'dFoF_first_sugar', 'dFoF_first_salt', 'intervals')
%save('trial_data_second.mat', 'dFoF_second','dFoF_second_inter', 'dFoF_second_sugar', 'dFoF_second_salt', 'intervals')

