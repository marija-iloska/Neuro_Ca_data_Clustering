function [store] = grouping(area_N, num_neurons)

neuron_idx = 1:num_neurons;
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
        store{ele} = [store{ele}, k1];

        % Change in area
        area_N(ele,:) = mean(area_N(store{ele}, :), 1);
    else
        % Store Group
        store{count} = [k1, k_new];
        group_k = [group_k, k1];

        % Change in area
        area_N(k1,:) = mean(area_N(store{count}, :), 1);
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

end