
function [train_data,test_data,class_2,class_1] = generate_sin_ts(ts_len, sin_len,n_samples)

    %default input n_sample = total number of data samples
    train_data = []; test_data = [];
    %     sin_len = 64;n_samples=4;
    x = (0:2*pi/sin_len:2*pi);
    sin_wave0 = 6*sin(x);
    sin_wave0 = sin_wave0(1:sin_len);
     %ts_c1 is the sythetic random walk + sine wave time series
    rnd_len = (ts_len-sin_len)/2;
    
    %generate train data
    for i=(1:n_samples/2)
        ts_c1 = [1];
        rnd = cumsum(randn(rnd_len,1));
        rnd_1 = rnd';
        ts_c1 = horzcat(ts_c1,rnd_1);
        %add Guassian noises to sine wave
        sin_wave = awgn(sin_wave0,10);
        %
        diff = ts_c1(end) - sin_wave(1);
        sin_wave = sin_wave + diff;
        ts_c1 = horzcat(ts_c1,sin_wave);
        rnd = cumsum(randn(rnd_len,1));
        rnd_1 = rnd';
        diff = ts_c1(end) - rnd_1(1);
        rnd_1 = rnd_1 + diff; % eliminate the gap between TS piece and random walk    
        ts_c1 = horzcat(ts_c1,rnd_1);

        train_data = vertcat(train_data,ts_c1);
        rnd_2 = cumsum(randn(ts_len,1));
        ts_c2 = [2,rnd_2'];
        train_data = vertcat(train_data,ts_c2);
    end
    
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
    
    %generate test data n_test = n_train/4
    n_test_c1 = mod(round(rand * 1000),300);
    n_test_c2 = 300 - n_test_c1;
    for i=1:n_test_c1
        ts_c1 = [1];
        rnd = cumsum(randn(rnd_len,1));
        rnd_1 = rnd';
        ts_c1 = horzcat(ts_c1,rnd_1);
        %add Guassian noises to sine wave
        sin_wave = awgn(sin_wave0,10);
        %
        diff = ts_c1(end) - sin_wave(1);
        sin_wave = sin_wave + diff;
        ts_c1 = horzcat(ts_c1,sin_wave);
        rnd = cumsum(randn(rnd_len,1));
        rnd_1 = rnd';
        diff = ts_c1(end) - rnd_1(1);
        rnd_1 = rnd_1 + diff; % eliminate the gap between TS piece and random walk    
        ts_c1 = horzcat(ts_c1,rnd_1);

        test_data = vertcat(test_data,ts_c1);

    end
    
    for i=1:n_test_c2
        rnd_2 = cumsum(randn(ts_len,1));
        ts_c2 = [2,rnd_2'];
        test_data = vertcat(test_data,ts_c2);
    end
    
    figure; subplot(211); plot(class_1(2:end)); title('Concatenated class 1'); 
    subplot(212); plot(class_2(2:end)); title('Concatenated class 2');
    disp(['Train size:', num2str(n_samples/2),' class A and ', num2str(n_samples/2),' class B']);
    disp(['Test size:',num2str(n_test_c1+n_test_c2),' in total']);
    disp('Finish preparing data');
%     
%     figure;
%     plot(train_data(1,:));
%     hold on;
%     plot(train_data(2,:));
end