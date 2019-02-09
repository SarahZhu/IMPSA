function [train_data, test_data, class_1, class_2, class_1_NoNaN, class_2_NoNaN] = prepare_data_tsv(dataset,n_train,n_test,prepare_mode)
    % prepare data for shapelet discovery
    % load train set, concatenate each class into a time series
    % add NaN at the concatenation point
    % datasets selected only have two classes
    % ***NOTE***: this is for train/test data with only 2 classes. Please
    % preprocess the data to 2-class dataset first
delimiterIn='\t';
train_data = importdata([dataset,'_TRAIN.tsv'],delimiterIn,0);

switch prepare_mode
    case 1
        disp('1.Using train data as many as possible with balanced classes');
        classes = unique(train_data(:, 1));
        class_1 = classes(1); class_1_NoNaN = classes(1);
        class_2 = classes(2); class_2_NoNaN = classes(2);
        count_cls1 = 0; count_cls2 = 0;
        [rows, ~] = size(train_data);
        for i = 1:rows
            if train_data(i, 1) == classes(1) 
                class_1 = [class_1, train_data(i, 2:end), NaN];
                class_1_NoNaN = [class_1_NoNaN,train_data(i, 2:end)];
                count_cls1 = count_cls1+1;
            elseif train_data(i, 1) == classes(2) 
                class_2 = [class_2, train_data(i, 2:end), NaN];
                class_2_NoNaN = [class_2_NoNaN,train_data(i, 2:end)];
                count_cls2 = count_cls2+1;
            end
        end
        
        %-----------balance #samples in the train set for each class-----------%

        len_cls1 = size(class_1); len_cls1_NoNaN = size(class_1_NoNaN);
        len_cls2 = size(class_2); len_cls2_NoNaN = size(class_2_NoNaN);

        if len_cls1(2) > len_cls2(2)
            class_1 = class_1(1:len_cls2(2));
            class_1_NoNaN = class_1_NoNaN(1:len_cls2_NoNaN(2));
        else
            class_2 = class_2(1:len_cls1(2));
            class_2_NoNaN = class_2_NoNaN(1:len_cls1_NoNaN(2));
        end

        %-----------balance #samples in the train set for each class-----------%
        disp(['Total #train samples: ',num2str(min(count_cls1,count_cls2)*2),' #train samples of class 1:',num2str(min(count_cls1,count_cls2)),' #train samples of class 2:',num2str(min(count_cls1,count_cls2))]);
    case 2
        disp('2.Using balanced train data for both two classes, the number of train data in total is specified by n_train');
        classes = unique(train_data(:, 1));
        class_1 = classes(1); class_1_NoNaN = classes(1);
        class_2 = classes(2); class_2_NoNaN = classes(2);
        count_cls1 = 0; count_cls2 = 0;
        [rows, ~] = size(train_data);
        for i = 1:rows
            if train_data(i, 1) == classes(1) && count_cls1 < n_train/2
                class_1 = [class_1, train_data(i, 2:end), NaN];
                class_1_NoNaN = [class_1_NoNaN,train_data(i, 2:end)];
                count_cls1 = count_cls1+1;
            elseif train_data(i, 1) == classes(2) && count_cls2 < n_train/2
                class_2 = [class_2, train_data(i, 2:end), NaN];
                class_2_NoNaN = [class_2_NoNaN,train_data(i, 2:end)];
                count_cls2 = count_cls2+1;
            end
        end
        
        if min(count_cls1,count_cls2) < n_train/2
            disp(['Insufficient train samples, the maximum of n_train is', num2str(min(count_cls1,count_cls2)),'. Please try again.']);
            exit;
        end
        
        disp(['Total #train samples: ',num2str(n_train),' #train samples of class 1:',num2str(count_cls1),' #train samples of class 2:',num2str(count_cls2)]);

    case 3
        disp('3.Using a selection of train dataset with a specified number(n_train) and balanced classes');
        train_data = train_data(1:n_train,:);

        classes = unique(train_data(:, 1));
        class_1 = classes(1); class_1_NoNaN = classes(1);
        class_2 = classes(2); class_2_NoNaN = classes(2);
        count_cls1 = 0; count_cls2 = 0;

        [rows, ~] = size(train_data);
        for i = 1:rows
            if train_data(i, 1) == classes(1)
                class_1 = [class_1, train_data(i, 2:end), NaN];
                class_1_NoNaN = [class_1_NoNaN,train_data(i, 2:end)]; %concatenating time series
                count_cls1 = count_cls1+1;
            elseif train_data(i, 1) == classes(2)
                class_2 = [class_2, train_data(i, 2:end), NaN];
                class_2_NoNaN = [class_2_NoNaN,train_data(i, 2:end)];
                count_cls2 = count_cls2+1;
            end
        end
        %-----------balance #samples in the train set for each class-----------%
        len_cls1 = size(class_1); len_cls1_NoNaN = size(class_1_NoNaN);
        len_cls2 = size(class_2); len_cls2_NoNaN = size(class_2_NoNaN);

        if len_cls1(2) > len_cls2(2)
            class_1 = class_1(1:len_cls2(2));
            class_1_NoNaN = class_1_NoNaN(1:len_cls2_NoNaN(2));
        else
            class_2 = class_2(1:len_cls1(2));
            class_2_NoNaN = class_2_NoNaN(1:len_cls1_NoNaN(2));
        end
        %-----------balance #samples in the train set for each class-----------%
        disp(['Total #train samples: ',num2str(min(count_cls1,count_cls2)*2),' #train samples of class 1:',num2str(min(count_cls1,count_cls2)),' #train samples of class 2:',num2str(min(count_cls1,count_cls2))]);      
    otherwise
        disp('Invalid mode for data preparation, please input 1,2 or 3 as prepare mode');
        exit;
end         
        

test_data = importdata([dataset,'_TEST.tsv'],delimiterIn,0);
test_data = test_data(1:n_test,:);

figure; subplot(211); plot(class_1(2:end)); title('Concatenated class 1'); 
subplot(212); plot(class_2(2:end)); title('Concatenated class 2');

disp('Finish preparing data');

end

% % Original Version
% train_data = train_data(1:n_train,:);
% 
% test_data = importdata([dataset,'_TEST.tsv'],delimiterIn,0);
% test_data = test_data(1:n_test,:);
% 
% classes = unique(train_data(:, 1));
% 
% class_1 = classes(1);
% class_2 = classes(2);
% 
% [rows, ~] = size(train_data);
% for i = 1:rows
%     if train_data(i, 1) == classes(1)
%         class_1 = [class_1, train_data(i, 2:end), NaN]; %concatenating time series
%     elseif train_data(i, 1) == classes(2)
%         class_2 = [class_2, train_data(i, 2:end), NaN];
%     end
% end
% 
% 
% % disp(['original class_1: ',num2str(size(class_1'))]);
% % disp(['original class_2: ',num2str(size(class_2'))]);
% 
% %-----------balance #samples in the train set for each class-----------%
% 
% len_cls1 = size(class_1);
% len_cls2 = size(class_2);
% 
% if len_cls1(2) > len_cls2(2)
%     class_1 = class_1(1:len_cls2(2));
% else
%     class_2 = class_2(1:len_cls1(2));
% end
% 
% %-----------balance #samples in the train set for each class-----------%
% 
% % disp(['balanced class_1: ',num2str(size(class_1))]);
% % disp(['balanced class_2: ',num2str(size(class_2))]);



