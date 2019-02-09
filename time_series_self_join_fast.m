% Compute the self-similarity join of a given time series A
% Chin-Chia Michael Yeh 01/26/2016
%
% [MatrixProfile, MatrixProfileIndex] = Time_series_Self_Join_Fast(A, SubsequenceLength)
% Output:
%     MatrixProfile: matrix profile of the self-join (vector)
%     MatrixProfileIndex: matrix profile index of the self-join (vector)
% Input:
%     A: input time series (vector)
%     SubsequenceLength: interested subsequence length (scalar)
%
% Chin-Chia Michael Yeh, Yan Zhu, Liudmila Ulanova, Nurjahan Begum, Yifei Ding, Hoang Anh Dau, Diego
% Furtado Silva, Abdullah Mueen, Eamonn Keogh. All Pairs Similarity Joins for Time Series Subsequences.
% SIGKDD 2016
%

function [MatrixProfile, MPindex] = time_series_self_join_fast(A, SubsequenceLength)
%% set trivial match exclusion zone
exclusionZone = round(SubsequenceLength/2);

%% check input
if SubsequenceLength > length(A)/2
    error('Error: Time series is too short relative to desired subsequence length');
end
if SubsequenceLength < 4
    error('Error: Subsequence length must be at least 4');
end
if length(A) == size(A, 2)
   A = A'; 
end

%% locate nan and inf
proLen = length(A) - SubsequenceLength + 1;
isSkip = false(proLen, 1);
for i = 1:proLen
    if any(isnan(A(i:i + SubsequenceLength - 1))) || ...
            any(isinf(A(i:i + SubsequenceLength - 1)))
        isSkip(i) = true;
    end
end

% replace NaN with 0s for the algoritm to work
A(isnan(A) | isinf(A)) = 0;

%% initialization
MatrixProfileLength = length(A) - SubsequenceLength + 1;
MatrixProfile = zeros(MatrixProfileLength, 1);
MPindex = zeros(MatrixProfileLength, 1);
[X, n, sumx2, sumx, meanx, sigmax2, sigmax] = ...
    fastfindNNPre(A, SubsequenceLength);

%% compute the matrix profile
pickedIdx = randperm(MatrixProfileLength);
for i = 1:MatrixProfileLength
    % compute the distance profile
    idx = pickedIdx(i);
    
    % skip subsequence containing NaN
    if isSkip(idx)
        continue
    end
    
    subsequence = A(idx:idx+SubsequenceLength-1);
    distanceProfile = fastfindNN(X, subsequence, n, SubsequenceLength, ...
        sumx2, sumx, meanx, sigmax2, sigmax);
    distanceProfile = abs(distanceProfile);
    
    % replace distance profile correspondings to subsequence contains NaN
    % with inf
    distanceProfile(isSkip) = inf;
    
    % apply exclusion zone
    exclusionZoneStart = max(1, idx-exclusionZone);
    exclusionZoneEnd = min(MatrixProfileLength, idx+exclusionZone);
    distanceProfile(exclusionZoneStart:exclusionZoneEnd) = inf;
    
    % figure out and store the nearest neighbor
    if i == 1
        MatrixProfile = distanceProfile;
        MPindex(:) = idx;
        [MatrixProfile(idx), MPindex(idx)] = min(distanceProfile);
    else
        updatePos = distanceProfile < MatrixProfile;
        MPindex(updatePos) = idx;
        MatrixProfile(updatePos) = distanceProfile(updatePos);
        [MatrixProfile(idx), MPindex(idx)] = min(distanceProfile);
    end
    
end

%% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
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