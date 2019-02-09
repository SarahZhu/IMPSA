function [train_data, test_data, class_2, class_1] = prepare_data(dataset,n_samples)
% prepare data for shapelet discovery
% load train set, concatenate each class into a time series
% add NaN at the concatenation point
% datasets selected only have two classes

train_data = load([dataset, '_train.mat']);
train_data = train_data.train_data;
% train_data = train_data.data;
test_data = load([dataset, '_test.mat']);
test_data = test_data.test_data;
% test_data = test_data.data;
% 
train_data = train_data(1:n_samples,:);
% test_data = test_data(1:4:n_samples,:);


classes = unique(train_data(:, 1));


class_1 = classes(1);
class_2 = classes(2);

[rows, ~] = size(train_data);
for i = 1:rows
    if train_data(i, 1) == classes(1)
        class_1 = [class_1, train_data(i, 2:end), NaN];
    elseif train_data(i, 1) == classes(2)
        class_2 = [class_2, train_data(i, 2:end), NaN];
    end
end

% disp(['original class_1: ',num2str(size(class_1'))]);
% disp(['original class_2: ',num2str(size(class_2'))]);

%-----------balance #samples in the train set for each class-----------%

len_cls1 = size(class_1);
len_cls2 = size(class_2);

if len_cls1(2) > len_cls2(2)
    class_1 = class_1(1:len_cls2(2));
else
    class_2 = class_2(1:len_cls1(2));
end

%-----------balance #samples in the train set for each class-----------%

% disp(['balanced class_1: ',num2str(size(class_1))]);
% disp(['balanced class_2: ',num2str(size(class_2))]);

figure; subplot(211); plot(class_1(2:end)); title('Concatenated class 1'); 
subplot(212); plot(class_2(2:end)); title('Concatenated class 2');

disp('Finish preparing data');

end

