clear all
close all
clc


% Choose session
session = {'s202', 's301', 's302', 's313', 'sf203', 'sf309', 'sf311'};

idx = 2;
ca_str = join(['Ca_Data/Ca_', session{idx}, '.mat']);

% Load data
load(ca_str)

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
%dFoF_norm = (dFoF -min_F)./w_norm;
% dFoF_norm = (dFoF - min_F);
% dFoF_norm = dFoF_norm./(max(dFoF_norm));
dFoF_norm = dFoF;


% JOIN all "taste" trials
for i = 1 : num_tastes

    % Which taste
    trial = trial_frames{i};
    range = [];
    first = [];
    second = [];
    clear psth_entire psth_first psth_second 

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
        psth_entire(n,1:num_neurons,:) = dFoF_norm(:,temp)- mean(dFoF_norm(:,temp(1:10)),2);
        psth_first(n, 1:num_neurons,:)  = dFoF_norm(:,first_temp)- mean(dFoF_norm(:,temp(1:10)),2);
        psth_second(n,1:num_neurons,:) =  dFoF_norm(:,second_temp) - mean(dFoF_norm(:, second_temp(1:10)),2);
    end

    % Get JOINT dFoF ranges
    entire_trial{i} = range;
    first_half{i} = first;
    second_half{i} = second;
   

    % PSTH
    psth_trial{i} = squeeze(mean( psth_entire, 1));
    psth_taste{i} = squeeze(mean( psth_first, 1));
    psth_decis{i} = squeeze(mean( psth_second, 1));


    % Name and get values using indices
    str = join(['dFoF_entire', '_', intervals{i}]);
    assignin('base',str, dFoF(:, entire_trial{i}));

    str = join(['dFoF_first', '_', intervals{i}]);
    assignin('base',str, dFoF(:, first_half{i}));

    str = join(['dFoF_second', '_', intervals{i}]);
    assignin('base',str, dFoF(:, second_half{i}));

end

for i = 1:num_tastes
    for n = 1:num_neurons
        psth_taste{i}(n,:) = psth_taste{i}(n,:) - min(psth_taste{i}(n,:));
        psth_taste{i}(n,:) = psth_taste{i}(n,:)/max(psth_taste{i}(n,:));

        psth_trial{i}(n,:) = psth_trial{i}(n,:) - min(psth_trial{i}(n,:));
        psth_trial{i}(n,:) = psth_trial{i}(n,:)/max(psth_trial{i}(n,:));

        psth_decis{i}(n,:) = psth_decis{i}(n,:) - min(psth_decis{i}(n,:));
        psth_decis{i}(n,:) = psth_decis{i}(n,:)/max(psth_decis{i}(n,:));
    end
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
filename = join(['Sort_Data/trial_', str_ses, '.mat']);
save(filename, 'dFoF_sugar', 'dFoF_salt', 'frames', 'split', ...
    'dFoF', 'psth_decis', 'psth_taste', 'psth_trial', 'intervals')


