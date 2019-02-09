function [] = find_sub_dist_for_each_class(train_data, shapelet)
% given a shapelet, find the subsequence distance for each class

% normalize shapelet
shapelet = zscore(shapelet);

classes = sort(unique(train_data(:, 1)));

class_1 = train_data(train_data(:, 1) == classes(1), :);
[rows_1, ~] = size(class_1);
class_2 = train_data(train_data(:, 1) == classes(2), :);
[rows_2, ~] = size(class_2);
class_1_sub_dists = zeros(rows_1, 1);
class_2_sub_dists = zeros(rows_2, 1);
for i = 1:rows_1
    ts = class_1(i, :);
    sub_dist = compute_subsequence_distance(shapelet, ts(2:end));
    class_1_sub_dists(i) = sub_dist;
end
for i = 1: rows_2
    ts = class_2(i, :);
    sub_dist = compute_subsequence_distance(shapelet, ts(2:end));
    class_2_sub_dists(i) = sub_dist;
end

class_1_sub_dists = sort(class_1_sub_dists);
class_2_sub_dists = sort(class_2_sub_dists);

class_1_time_vector = ones(rows_1, 1);
class_2_time_vector = 2 * ones(rows_2, 1);

figure; plot(class_1_time_vector, class_1_sub_dists, 'x',...
    'color', 'r', 'LineWidth', 1);
hold on; plot(class_2_time_vector, class_2_sub_dists, 'o', ...
    'color', 'b', 'LineWidth', 1);
title('Sorted subsequence distance for each class');
xlim([0 4]);
ylim([0 7]);
legend('Class 1', 'Class 2');

end



   