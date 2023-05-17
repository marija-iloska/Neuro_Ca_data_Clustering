clear all
close all
clc

% Load data
load Ca_data.mat
load trial_data.mat

% SUGAR, SALT, INTER
% Compute absolute area for each Neuron
for n = 1:num_neurons
    area_sugar(n) = trapz(abs(dFoF_sugar(n,:)));
    area_salt(n) =  trapz(abs(dFoF_salt(n,:)));
    area_inter(n) =  trapz(abs(dFoF_inter(n,:)));

    % Sort area and get indices
    [areaSU, idx_SU] = sort(area_sugar, 'descend');
    [areaSA, idx_SA] = sort(area_salt, 'descend');
    [areaIN, idx_IN] = sort(area_inter, 'descend');
end

