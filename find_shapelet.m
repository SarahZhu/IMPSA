function [shapelet_idxs, shapelet_label] = find_shapelet(A, B, sub_len)
% find the shapelet by the Matrix Profile approach
% A: the concatenation of class 1
% B: the concatenation of class 2

number_of_shapelets = 2;
shapelet_idxs = zeros(number_of_shapelets, 1);

shapelet_label = A(1);

A = A(2:end);
B = B(2:end);

[MPAA, ~] = time_series_self_join_fast(A, sub_len);
[MPAB, ~] = time_series_join_fast(A, B, sub_len);

MP_diff = MPAB - MPAA;
MP_diff_4_plot = MP_diff;
MP_diff(~isfinite(MP_diff))=0;
MP_diff(isnan(MP_diff)) = 0;
[sorted_MP_diff, sort_idx] = sort(MP_diff, 'descend');
sort_idx(isnan(sorted_MP_diff)) = [];

for i = 1: number_of_shapelets
    shapelet_idxs(i) = sort_idx(i);
end

% % plot the most potential shapelet so far
% [~, shapelet_index] = max(MP_diff);
% shapelet = A(shapelet_index:shapelet_index + sub_len - 1);
% idx_of_ts = round(shapelet_index / 151);
% ts = A(151 * (idx_of_ts - 1) + 1: 151 * idx_of_ts);
% shapelet_local_idx = mod(shapelet_index, 151);
% shapelet_time_vector = shapelet_local_idx:shapelet_local_idx + sub_len - 1;
% figure; plot(ts); hold on;  
% plot(shapelet_time_vector, shapelet, 'color', rand(3,1), 'LineWidth', 3);
% title('Time series and shapelet part');
% 
% figure; plot(shapelet, 'color', rand(3, 1)); title('Best shapelet');

shapelet_idxs_values = MP_diff(shapelet_idxs);
%----------------------------------------------------------------------------------
figure; plot(MP_diff_4_plot, 'color', rand(3, 1)); title('Matrix Profile difference');
hold on; plot(shapelet_idxs, shapelet_idxs_values, 'o');

figure; plot(MPAA, 'color', rand(3, 1), 'LineWidth', 2); 
hold on; plot(MPAB, 'color', rand(3, 1)); title('Matrix Profile');
%----------------------------------------------------------------------------------
disp('Finish finding shapelet');

end