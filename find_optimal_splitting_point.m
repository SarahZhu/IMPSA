function [optimal_splitting_point, bsf_accuracy] = ...
    find_optimal_splitting_point(train_data, shapelet, shapelet_label,dist_mode)
% given a shapelet, find the optimal splitting point

% normalize shapelet
if strcmp(dist_mode,'Original') == 1
    shapelet = zscore(shapelet);
end
    

[rows, ~] = size(train_data);
sub_dists = zeros(rows, 1);

for i = 1:rows
    ts = train_data(i, 2:end);
    sub_dists(i) = compute_subsequence_distance(shapelet, ts, dist_mode);
    
end

[sub_dists, sorted_idx] = sort(sub_dists);
true_label = train_data(:, 1);
true_label = true_label(sorted_idx);

potential_split_points = movsum(sub_dists, 2) ./ 2;

% M = movsum(___,dim) returns the array of 
% sliding sums along dimension dim for
% any of the previous syntaxes. 
% For example, if A is a matrix, 
%     then movsum(A,k,2) operates 
%     along the columns of A, computing the k-element sliding sum for each row.

% each subsequence distance is a potential best spliting point
% test them one by one
tol = eps(0.5);
bsf_accuracy = 0;
optimal_splitting_point = NaN;

for i = 1:length(potential_split_points)
   count = 0;
   candidate = potential_split_points(i);
   
   for j = 1:rows
       if sub_dists(j) <= candidate 
           if true_label(j) == shapelet_label
           count = count + 1;
           end
       else
           if true_label(j) ~= shapelet_label
               count = count + 1;
           end
       end
   end
   
   train_accuracy = count / rows;
  
   if train_accuracy > bsf_accuracy
       bsf_accuracy = train_accuracy;
       optimal_splitting_point = candidate;
   end
   if abs(bsf_accuracy-1.0)< tol  %Added by Sarah to save time
       break
   end  
end

fprintf('Optimal splitting point = %f \n', optimal_splitting_point);
fprintf('Best train accuracy = %f \n', bsf_accuracy);
end



   