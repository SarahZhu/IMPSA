function [bsf_dist] = compute_subsequence_distance(subsequence, time_series, mode)
% compute the subsequence distance between a subsequence and a time series
% assume that the subsequence has already been normalized

    if length(subsequence) > length(time_series)
       error('Time series must be longer than subsequence length'); 
    end
    switch mode
        case 'Original'
            bsf_dist = inf;
            subsequence_count = length(time_series) - length(subsequence) + 1;
            for i = 1: subsequence_count
                ts_sub = time_series(i: i + length(subsequence) - 1);
                ts_sub = zscore(ts_sub);
                sub_dist = pdist2(subsequence, ts_sub);
                if bsf_dist >  sub_dist
                    bsf_dist = sub_dist;
                end
            end 
        case 'MASS'
            [dist_profile] = MASS_V2(time_series,subsequence);
            [sorted_dist_profile,~] = sort(dist_profile, 'ascend');
            bsf_dist = sorted_dist_profile(1);
        otherwise
            disp('Invalid mode, please try again (mode can be set as Original or MASS)');
            exit;
end