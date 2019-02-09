function [shapelet_idxs, shapelet_label] = find_shapelet2(A, B, sub_len)
% find the shapelet by the Matrix Profile approach
% A: the concatenation of class 1
% B: the concatenation of class 2

number_of_shapelets = 10;
shapelet_idxs = zeros(number_of_shapelets, 1);

shapelet_label = A(1);

A = A(2:end);
B = B(2:end);

% [MPAA, ~] = time_series_self_join_fast(A, sub_len);
[MPAA, ~] = PreSCRIMP(A, sub_len, 0.25);
% [MPAB, ~] = time_series_join_fast(A, B, sub_len);
[MPAB, ~] = PreSCRIMP_joinAB(A, B, sub_len, 0.25);
% global MP_diff sorted_MP_diff sort_idx
MP_diff = MPAB - MPAA;
% MP_diff_4_plot = MP_diff;
MP_diff(~isfinite(MP_diff))=0;
MP_diff(isnan(MP_diff)) = 0;
[sorted_MP_diff, sort_idx] = sort(MP_diff, 'descend');
sort_idx(isnan(sorted_MP_diff)) = [];

for i = 1: number_of_shapelets
    shapelet_idxs(i) = sort_idx(i);
end

disp('Finish finding shapelet');

end