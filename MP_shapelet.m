function [] = MP_shapelet()
dataset_name = 'data/GunPoint';
sub_len = 38;

[train_data, test_data, A, B] = prepare_data(dataset_name);

[shapelet_idxs, shapelet_label] = ...
    find_shapelet(B, A, sub_len);

% test each shapelet
B = B(2:end);
bsf_shapelet = [];
bsf_threshold = NaN;
bsf_accuracy = 0;
bsf_shapelet_idx = NaN;
train_accuracies = zeros(1, length(shapelet_idxs));
test_accuracies = zeros(1, length(shapelet_idxs));
for i = 1:length(shapelet_idxs)
    idx = shapelet_idxs(i);
    shapelet = B(idx:idx + sub_len - 1);
    [optimal_splitting_point, best_accuracy] = ...
    find_optimal_splitting_point(train_data, shapelet, shapelet_label);
    if bsf_accuracy < best_accuracy
        bsf_accuracy = best_accuracy;
        bsf_threshold = optimal_splitting_point;
        bsf_shapelet = shapelet;
        bsf_shapelet_idx = idx;
    end
    train_accuracies(i) = best_accuracy;
    test_accuracies(i) = run_shapelet_classification(test_data, ...
        shapelet, shapelet_label, optimal_splitting_point);
end

% test the best shapelet
run_shapelet_classification(test_data, bsf_shapelet, shapelet_label, ...
    bsf_threshold);

plot_best_shapelet(B, bsf_shapelet_idx);

figure; plot(train_accuracies, 'color', rand(3,1));
hold on; plot(test_accuracies, 'color', rand(3, 1), 'LineWidth', 3);

% plot the subsequence distance for train data
find_sub_dist_for_each_class(train_data, bsf_shapelet);
% plot the subsequence distance for test data
find_sub_dist_for_each_class(test_data, bsf_shapelet);

end