function [] = plot_best_shapelet2(A,shapelet_index,sub_len,ts_len)
% plot the time series and shapelet
% sub_len = 38;
% ts_len = 151;
% A = A(2:end);
shapelet = A(shapelet_index:shapelet_index + sub_len - 1);
%idx_of_ts = round(shapelet_index / ts_len);
idx_of_ts = ceil(shapelet_index / ts_len);
ts = A(ts_len * (idx_of_ts - 1) + 1: ts_len * idx_of_ts);
shapelet_local_idx = mod(shapelet_index, ts_len);
shapelet_time_vector = shapelet_local_idx:shapelet_local_idx + sub_len - 1;
figure; plot(ts); hold on;  
plot(shapelet_time_vector, shapelet, 'color', rand(3,1), 'LineWidth', 3);
title('Time series and shapelet part');

figure; plot(shapelet, 'color', rand(3, 1)); title('Best shapelet');
end