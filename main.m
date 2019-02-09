%% test MP shapelet
% % Arrowhead
% sub_len = 95;              % subsequence length
% ts_len = 251;

% %Gun
% sub_len = 38;              % subsequence length
% ts_len = 150;

% % Yoga
% sub_len = 88;              % subsequence length
% ts_len = 426;

% %Mallat
% sub_len = 635;              % subsequence length 530, 635
% ts_len = 1024;

% % StarLightCurves
sub_len = 326;              % subsequence length
ts_len = 1024;

%----------------------------------------------------------------------------------
% dataset_name = 'data/Gun'; % folder contains train and test data(mat format)
% [train_data, test_data, A, B] = prepare_data(dataset_name);
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
dataset_name = 'data/tsv/StarLightCurves';% folder contains train and test data(tsv format)
n_train = 10;
n_test = 300;
prepare_mode = 2;
[train_data, test_data, A, B] = prepare_data_tsv(dataset_name,n_train,n_test,prepare_mode);

%----------------------------------------------------------------------------------
% % Adding white gaussian noises, to and to avoid the NaN problem increase accuracy
% train_data = awgn(train_data,10);
% test_data = awgn(test_data,10);
%----------------------------------------------------------------------------------
%-------------------------------prepare mode---------------------------------------
% '1.Using train data as many as possible with balanced classes'
% '2.Using balanced train data for both two classes, the number of train data for each class is specified by n_train'
% '3.Using a selection of train dataset with a specified number(n_train) and balanced classes'
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------

[shapelet_idxs, shapelet_label] = ...
    find_shapelet(B, A, sub_len);


% [shapelet_idxs, shapelet_label] = ...
%     find_shapelet(A, B, sub_len);

%test each shapelet and update best-so-far shapelet
B = B(2:end);
% A = A(2:end);
bsf_shapelet = [];
bsf_threshold = NaN;
bsf_accuracy = 0;
bsf_shapelet_idx = NaN;
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
end

run_shapelet_classification(test_data, bsf_shapelet, shapelet_label, ...
    bsf_threshold);
% 
% plotMPnMPDiff(B, A, sub_len);

%plot_best_shapelet(B, bsf_shapelet_idx)
% plot_best_shapelet2(A, bsf_shapelet_idx,sub_len,ts_len);
plot_best_shapelet2(B, bsf_shapelet_idx,sub_len,ts_len)
% % plot the subsequence distance vector for each class in training
% find_sub_dist_for_each_class(train_data, shapelet);
% 
% % plot the subsequence distance vector for each class in testing
% find_sub_dist_for_each_class(test_data, shapelet);
% 
% % find optimal spliting point on train data
% [optimal_splitting_point, best_accuracy] = ...
%     find_optimal_splitting_point(train_data, shapelet, shapelet_label);
% 
% % compute classification accuracy on test data
% run_shapelet_classification(test_data, shapelet, ...
%     shapelet_label, optimal_splitting_point);
% 
% %% test lexiang shapelet
% shapelet_label = 1;
% idx = 8; % not 45 as the author said
% start_pos = 101;
% ts = train_data(idx, 2:end);
% time_vec = start_pos:start_pos+sub_len-1;
% lexiang_shapelet = ts(start_pos:start_pos+sub_len-1);
% 
% lexiang_threshold = 5.7226; % not sqrt(38.94)
% run_shapelet_classification(test_data, lexiang_shapelet, ...
%     shapelet_label, lexiang_threshold);

% figure; plot(ts); hold on; plot(time_vec, lexiang_shapelet, 'LineWidth', 3);