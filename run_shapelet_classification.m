function [accuracy] = run_shapelet_classification(test_data, ...
    shapelet, shapelet_label, optimal_splitting_point,mode)

% normalize shapelet
shapelet = zscore(shapelet);

[rows, ~] = size(test_data);
true_label = test_data(:, 1);
count = 0;
for i = 1:rows
   ts = test_data(i, 2:end);
   sub_dist = compute_subsequence_distance(shapelet, ts,mode);
   if sub_dist <= optimal_splitting_point
       if true_label(i) == shapelet_label
           count = count + 1;
       end
   else
      if true_label(i) ~= shapelet_label
          count = count + 1;
      end
   end
end

accuracy = count / rows;
fprintf('Test accuracy = %f \n', accuracy);
end