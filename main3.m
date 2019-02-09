%% This program is for plotting training and testing accuracy

%% test MP shapelet

% % Arrowhead 
% % sub_len = 44;               % subsequence length
% ts_len = 251;

% %Gun
% sub_len = 38;                 % subsequence length
% ts_len = 150;

% % Yoga
% % sub_len = 200;              % subsequence length
% ts_len = 426;

% %Mallat
% % sub_len = 450;              % subsequence length
% ts_len = 1024;

% %StarLightCurves
% % sub_len = 300;              % subsequence length
% ts_len = 1024;


dist_mode='MASS';
% dist_mode='Original';% method used to compute subsequence distance = {'Original','MASS'}
%----------------------------------------------------------------------------------
% dataset_name = 'data/FuzzySine128'; % folder contains train and test data(mat format)
%  n_train = 128; ts_len = 1024; sl_temp = 100;
% [train_data, test_data, A, B] = prepare_data(dataset_name,n_train);
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------

% dataset_name = 'data/tsv/StarLightCurves_c1c2';% folder contains train and test data(tsv format)
% n_train = 128; ts_len = 1024; sl_temp = 300;
% n_test = 300;
% prepare_mode = 2;
% [train_data, test_data, A, B,~,~] = prepare_data_tsv(dataset_name,n_train,n_test,prepare_mode);
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
dataset_name = 'data/tsv/TwoPatterns';% folder contains train and test data(tsv format)
n_train = 250; ts_len = 128; %sl_temp = 38;
n_test = 1000;
prepare_mode = 2;
[train_data, test_data, A, B,~,~] = prepare_data_tsv(dataset_name,n_train,n_test,prepare_mode);
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------

% dataset_name = 'data/tsv/GunPoint';% folder contains train and test data(tsv format)
% n_train = 32; ts_len = 150; sl_temp = 38;
% n_test = 150;
% prepare_mode = 2;
% [train_data, test_data, A, B,~,~] = prepare_data_tsv(dataset_name,n_train,n_test,prepare_mode);
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------

% dataset_name = 'data/tsv/Mallat_c3c4';% folder contains train and test data(tsv format)
% n_train = 100; ts_len = 1024; sl_temp = 300;
% n_test = 300;
% prepare_mode = 2;
% [train_data, test_data, A, B] = prepare_data_tsv(dataset_name,n_train,n_test,prepare_mode);
%----------------------------------------------------------------------------------
% sin_len = 64;
% [train_data,test_data,A,B] = generate_sin_ts(ts_len, sin_len,n_samples);
%----------------------------------------------------------------------------------
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
disp('Start training, profiler on');
profile on;
tol = eps(0.5);

rec_acc_test = []; rec_acc_train = [];
bsf_glb_accuracy = 0;
bsf_glb_shapelet = [];
bsf_glb_threshold = NaN;
bsf_glb_shapelet_idx = NaN;
sub_len_glb = 0;
step = round(0.05*ts_len);
for sub_len = 4:step:ts_len %sl_temp:1:sl_temp       %sub_len = 4:step:ts_len
    disp(['The current sub_len is ',num2str(sub_len),':']);
    %----------------------------------------------------------------------------------
    [shapelet_idxs, shapelet_label] = ...
        find_shapelet2(B, A, sub_len);


    % [shapelet_idxs, shapelet_label] = ...
    %     find_shapelet(A, B, sub_len);

    %test each shapelet and update best-so-far shapelet  n
    B_p = B(2:end);
    % A = A(2:end);
    bsf_shapelet = [];
    bsf_threshold = NaN;
    bsf_accuracy = 0;
    bsf_shapelet_idx = NaN;
    for i = 1:length(shapelet_idxs)
        idx = shapelet_idxs(i);
        shapelet = B_p(idx:idx + sub_len - 1);
        [optimal_splitting_point, best_accuracy] = ...
        find_optimal_splitting_point(train_data, shapelet, shapelet_label, dist_mode);    

        if bsf_accuracy < best_accuracy
            bsf_accuracy = best_accuracy;
            bsf_threshold = optimal_splitting_point;
            bsf_shapelet = shapelet;
            bsf_shapelet_idx = idx;
        end
%         if abs(bsf_accuracy-1.0)< tol %Added by Sarah to save
%         time,comment out when plotting train&test accuracy
%            break
%         end
    end
   
    rec_acc_train = [rec_acc_train,bsf_accuracy];
    
    acc_test = run_shapelet_classification(test_data, bsf_shapelet, shapelet_label, ...
    bsf_threshold, dist_mode);
    rec_acc_test = [rec_acc_test, acc_test]; 
    
    
    %update global best shapelet info
    if bsf_glb_accuracy < bsf_accuracy %changed from < to <= here
        bsf_glb_accuracy = bsf_accuracy;
        bsf_glb_shapelet = bsf_shapelet;
        bsf_glb_threshold = bsf_threshold;
        bsf_glb_shapelet_idx = bsf_shapelet_idx;
        sub_len_glb = sub_len;
    end
%     if abs(bsf_glb_accuracy-1.0)< tol %Added by Sarah to save time,comment out when plotting train&test accuracy
%        break
%     end
end

disp('Plotting test & train accuracy vs. sublength of shapelet candidates...');
figure;
x = 4:step:ts_len;
disp([length(x),length(rec_acc_train),length(rec_acc_test)]);
plot(x,rec_acc_train,'-o');hold on;plot(x,rec_acc_test,'-x');
legend({'Train','Test'},'Location','southwest');
title('Train and Test Accuracy vs. Lengths of Shapelet Candidates');
xlabel('Sublength') ;
ylabel('Accuracy') ;
for i=1:length(x)
    text(x(i)+3,rec_acc_train(i),num2str(rec_acc_train(i)),'Color','b','FontSize',14);
    text(x(i)-3,rec_acc_test(i),num2str(rec_acc_test(i)),'Color','r','FontSize',14);
end
hold off;

disp('Search of best sublength of shapelets finished');
disp(['The best sub_len is ',num2str(sub_len_glb)]);

profile viewer;
profsave;

disp(['Total number of training samples: ',num2str(n_train)]);
disp('Start classifying test data, profiler saved');
run_shapelet_classification(test_data, bsf_glb_shapelet, shapelet_label, ...
    bsf_glb_threshold, dist_mode);
disp('Plotting results...');
plotMPnMPDiff(B, A, sub_len_glb);

plot_best_shapelet2(B_p, bsf_glb_shapelet_idx,sub_len_glb,ts_len+1)

