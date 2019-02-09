function [PreMatrixProfile, PreMPindex] = PreSCRIMP(A, SubsequenceLength, step)
% Compute the approximate self similarity join of time series A with
% PreSCRIMP
% Matrix Pro?le XI: SCRIMP++: Time Series Motif Discovery at Interactive Speed. Yan Zhu, Chin-Chia Michael Yeh, Zachary Zimmerman, Kaveh Kamgar and Eamonn Keogh, ICDM 2018. 
% https://www.cs.ucr.edu/~eamonn/SCRIMP_ICDM_camera_ready_updated.pdf
% For details of the SCRIMP++ algorithm, see:
% "SCRIMP++: Motif Discovery at Interactive Speeds", submitted to ICDM 2018.
% Usage:
% [PreMatrixProfile, PreMPindex] = PreSCRIMP(A, SubsequenceLength, step)
% Output:
%     PreMatrixProfile: running matrix profile after PreSCRIMP is completed
%     (vector)
%     PreMPindex: running matrix profile index after PreSCRIMP is completed
%     (vector)
% Input:
%     A: input time series (vector)
%     SubsequenceLength: interested subsequence length (scalar)
%     step: s/SubsequenceLength, the step size of PreSCRIMP. For all
%     experiments in the paper, this input is set to a fixed value 0.25.

%% set trivial match exclusion zone
exclusionZone = round(SubsequenceLength/4);

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
%-------------------added by Sarah, handle NaN for shapelet algorithm----------------------------------

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
current_step = floor(SubsequenceLength*step);
pickedIdx = 1:current_step:MatrixProfileLength;

dotproduct = zeros(MatrixProfileLength,1);
refine_distance = inf(MatrixProfileLength,1);
orig_index = 1:MatrixProfileLength;
m = SubsequenceLength;

disp('Running PreSCRIMP-time series self join P_AA:');

for i = 1:length(pickedIdx)
    % compute the distance profile
    idx = pickedIdx(i);
    subsequence = A(idx:idx+SubsequenceLength-1);
    [distanceProfile, ~, ~, ~, ~] = fastfindNN(X, subsequence, n, SubsequenceLength, ...
        sumx2, sumx, meanx, sigmax2, sigmax);
    distanceProfile = abs(distanceProfile);
    
    %-------------------added by Sarah----------------------------------
    % replace distance profile correspondings to subsequence contains NaN
    % with inf
    distanceProfile(isSkip) = inf; 
    %-------------------------------------------------------------------
    
    % apply exclusion zone
    exclusionZoneStart = max(1, idx-exclusionZone);
    exclusionZoneEnd = min(MatrixProfileLength, idx+exclusionZone);
    distanceProfile(exclusionZoneStart:exclusionZoneEnd) = inf;
    
    % figure out and store the neareest neighbor
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
    
    idx_nn = MPindex(idx);
    idx_diff = idx_nn-idx;
    dotproduct(idx) = (m - MatrixProfile(idx)^2/2)*sigmax(idx)*sigmax(idx_nn) + m*meanx(idx)*meanx(idx_nn);
    
    endidx = min([MatrixProfileLength, idx+current_step-1, MatrixProfileLength-idx_diff]);
    dotproduct(idx+1:endidx) =  dotproduct(idx)+cumsum(A(idx+m:endidx+m-1).*A(idx_nn+m:endidx+m-1+idx_diff)-A(idx:endidx-1).*A(idx_nn:endidx-1+idx_diff));
    refine_distance(idx+1:endidx) = sqrt(abs(2*(m-(dotproduct(idx+1:endidx)-m*meanx(idx+1:endidx).*meanx(idx_nn+1:endidx+idx_diff))./(sigmax(idx+1:endidx).*sigmax(idx_nn+1:endidx+idx_diff)))));
    
    beginidx = max([1, idx-current_step+1, 1-idx_diff]);
    dotproduct(idx-1:-1:beginidx) = dotproduct(idx)+cumsum(A(idx-1:-1:beginidx).*A(idx_nn-1:-1:beginidx+idx_diff)-A(idx-1+m:-1:beginidx+m).*A(idx_nn-1+m:-1:beginidx+idx_diff+m));
    refine_distance(beginidx:idx-1) = sqrt(abs(2*(m-(dotproduct(beginidx:idx-1)-m*meanx(beginidx:idx-1).*meanx(beginidx+idx_diff:idx_nn-1))./(sigmax(beginidx:idx-1).*sigmax(beginidx+idx_diff:idx_nn-1)))));
    
    update_pos1 = find(refine_distance(beginidx:endidx)<MatrixProfile(beginidx:endidx));
    MatrixProfile(update_pos1+beginidx-1) = refine_distance(update_pos1+beginidx-1);
    MPindex(update_pos1+beginidx-1) = orig_index(update_pos1+beginidx-1)+idx_diff;
    
    update_pos2 = find(refine_distance(beginidx:endidx)<MatrixProfile(beginidx+idx_diff:endidx+idx_diff));
    MatrixProfile(update_pos2+beginidx+idx_diff-1) = refine_distance(update_pos2+beginidx-1);
    MPindex(update_pos2+beginidx+idx_diff-1) = orig_index(update_pos2+beginidx+idx_diff-1)-idx_diff;
    
    if mod(i,1000)==0
        completed_perc = i/length(pickedIdx)*100;
        disp(['progress:', num2str(completed_perc), '%']);
    end
end

PreMatrixProfile = MatrixProfile;
PreMPindex = MPindex;
disp('PreSCRIMP-time series self join finished');

% m is winSize
function [X, n, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(x, m)
n = length(x);
%x(n+1:2*n) = 0;
X = fft(x);
cum_sumx = cumsum(x);
cum_sumx2 =  cumsum(x.^2);
sumx2 = cum_sumx2(m:n)-[0;cum_sumx2(1:n-m)];
sumx = cum_sumx(m:n)-[0;cum_sumx(1:n-m)];
meanx = sumx./m;
sigmax2 = (sumx2./m)-(meanx.^2);
sigmax = sqrt(sigmax2);

% m is winSieze
function [dist lastz dropval, sumy sumy2] = fastfindNN(X, y, n, m, sumx2, sumx, meanx, sigmax2, sigmax)
%x is the data, y is the query
%y = (y-mean(y))./std(y,1);                      %Normalize the query
dropval=y(1);
y = y(end:-1:1);                                %Reverse the query
y(m+1:n) = 0;

%The main trick of getting dot products in O(n log n) time
Y = fft(y);
Z = X.*Y;
z = ifft(Z);

%compute y stats -- O(n)
sumy = sum(y);
sumy2 = sum(y.^2);
meany=sumy/m;
sigmay2 = sumy2/m-meany^2;
sigmay = sqrt(sigmay2);

%computing the distances -- O(n) time
%dist = (sumx2 - 2*sumx.*meanx + m*(meanx.^2))./sigmax2 - 2*(z(m:n) - sumy.*meanx)./sigmax + sumy2;
%dist = 1-dist./(2*m);

dist = 2*(m-(z(m:n)-m*meanx*meany)./(sigmax*sigmay));
dist = sqrt(dist);
lastz=real(z(m:n));