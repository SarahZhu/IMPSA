% Compute the similarity join of a given time series A with another time series B
% Chin-Chia Michael Yeh 01/26/2016
%
% [MatrixProfile, MatrixProfileIndex] = Time_series_Join_Fast(A, B, SubsequenceLength) 
% Output:
%     MatrixProfile: matrix porfile of the join (vector)
%     MatrixProfileIndex: matrix porfile index of the join (vector)
% Input:
%     A: first input time series (vector)
%     B: second input time series (vector)
%     SubsequenceLength: interested subsequence length (scalar)
%
% Chin-Chia Michael Yeh, Yan Zhu, Liudmila Ulanova, Nurjahan Begum, Yifei Ding, Hoang Anh Dau, Diego
% Furtado Silva, Abdullah Mueen, Eamonn Keogh. All Pairs Similarity Joins for Time Series Subsequences.
% SIGKDD 2016
%

function [MatrixProfile, MPindex] = time_series_join_fast(A, B, SubsequenceLength) 
%% check input
if SubsequenceLength > length(A)/2 
    error('Error: Time series A is too short relative to desired subsequence length');
end
if SubsequenceLength > length(B)/2 
    error('Error: Time series B is too short relative to desired subsequence length');
end
if SubsequenceLength < 4 
    error('Error: Subsequence length must be at least 4');
end
if length(A) == size(A, 2)
   A = A'; 
end
if length(B) == size(B, 2)
   B = B'; 
end

%% locate nan and inf
ASubsequenceCount = length(A) - SubsequenceLength + 1;
BSubsequenceCount = length(B) - SubsequenceLength + 1;
isSkipA = false(ASubsequenceCount, 1);
isSkipB = false(BSubsequenceCount, 1);
for i = 1:ASubsequenceCount
    if any(isnan(A(i:i + SubsequenceLength - 1))) || ...
            any(isinf(A(i:i + SubsequenceLength - 1)))
        isSkipA(i) = true;
    end
end

for i = 1:BSubsequenceCount
    if any(isnan(B(i:i + SubsequenceLength - 1))) || ...
            any(isinf(B(i:i + SubsequenceLength - 1)))
        isSkipB(i) = true;
    end
end
% replace NaN with 0s for the algoritm to work
A(isnan(A) | isinf(A)) = 0;
B(isnan(B) | isinf(B)) = 0;

%% initialization
% ASubsequenceCount = length(A) - SubsequenceLength + 1;
% BSubsequenceCount = length(B) - SubsequenceLength + 1;
MatrixProfile = inf(ASubsequenceCount, 1);
MPindex = zeros(ASubsequenceCount, 1); 
[X, n, sumx2, sumx, meanx, sigmax2, sigmax] = ...
    fastfindNNPre(A, SubsequenceLength);

%% compute the matrix profile
pickedIdx = randperm(BSubsequenceCount);
for i = 1:BSubsequenceCount
    % compute the distance profile
    idx = pickedIdx(i);
    
    if isSkipB(idx)
        continue
    end
    
    subsequence = B(idx:idx+SubsequenceLength-1);
    distanceProfile = fastfindNN(X, subsequence, n, SubsequenceLength, ...
        sumx2, sumx, meanx, sigmax2, sigmax);
    distanceProfile = abs(distanceProfile);
    
    distanceProfile(isSkipA) = inf;
	
	% figure out and store the nearest neighbor
    if i == 1
        MatrixProfile = distanceProfile;
        MPindex(:) = idx;
    else
        updatePos = distanceProfile <= MatrixProfile;
        MPindex(updatePos) = idx;
        MatrixProfile(updatePos) = distanceProfile(updatePos);
    end
    
end

%% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
%  https://www.cs.unm.edu/~mueen/findNN.html(MASS)
function [X, n, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(x, m)
n = length(x);
x(n+1:2*n) = 0;
X = fft(x);
cum_sumx = cumsum(x);
cum_sumx2 =  cumsum(x.^2);
sumx2 = cum_sumx2(m:n)-[0;cum_sumx2(1:n-m)];
sumx = cum_sumx(m:n)-[0;cum_sumx(1:n-m)];
meanx = sumx./m;
sigmax2 = (sumx2./m)-(meanx.^2);
sigmax = sqrt(sigmax2);

function dist = fastfindNN(X, y, n, m, sumx2, sumx, meanx, sigmax2, sigmax)
%x is the data, y is the query
y = (y-mean(y))./std(y,1);                      %Normalize the query
y = y(end:-1:1);                                %Reverse the query
y(m+1:2*n) = 0;

%The main trick of getting dot products in O(n log n) time
Y = fft(y);
Z = X.*Y;
z = ifft(Z);

%compute y stats -- O(n)
sumy = sum(y);
sumy2 = sum(y.^2);

%computing the distances -- O(n) time
dist = (sumx2 - 2*sumx.*meanx + m*(meanx.^2))./sigmax2 - 2*(z(m:n) - sumy.*meanx)./sigmax + sumy2;
%dist = 1-dist./(2*m);
dist = sqrt(dist);