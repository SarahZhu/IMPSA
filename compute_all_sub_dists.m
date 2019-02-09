function [sub_dists] = compute_all_sub_dists (train_data, test_data, shapelet)

shapelet = zscore(shapelet);

[train_rows, ~] = size(train_data);
[test_rows, ~] = size(test_data);

sub_dists = zeros(train_rows + test_rows, 1);
for i = 1:train_rows
   ts = train_data(i, 2:end);
   sub_dist = compute_subsequence_distance(shapelet, ts);
   sub_dists(i) = sub_dist;
end

for i = 1:test_rows
   ts = test_data(i, 2:end);
   sub_dist = compute_subsequence_distance(shapelet, ts);
   sub_dists(i + train_rows) = sub_dist;
end

figure; plot(sort(sub_dists)); title('Sorted subsequence distance');

end