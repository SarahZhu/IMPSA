function [PreMatrixProfile, PreMPindex] = PreSCRIMP_joinAB(A, B, SubsequenceLength, step)
% Compute the approximate self similarity join of time series A with
% PreSCRIMP
% Author information omitted for ICDM review.
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


%% check input
if SubsequenceLength > length(A)/2
    error('Error: Time series is too short relative to desired subsequence length');
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
%-------------------added by Sarah----------------------------------


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
% MatrixProfileLength = length(A) - SubsequenceLength + 1;
% MatrixProfile = zeros(MatrixProfileLength, 1);
% MPindex = zeros(MatrixProfileLength, 1);
% [X, n, sumx2, sumx, meanx, sigmax2, sigmax] = ...
%     fastfindNNPre(A, SubsequenceLength);

% ASubsequenceCount = length(A) - SubsequenceLength + 1;
% BSubsequenceCount = length(B) - SubsequenceLength + 1;
MatrixProfile = inf(ASubsequenceCount, 1);
MPindex = zeros(ASubsequenceCount, 1); 
[X, n, sumx2, sumx, meanx, sigmax2, sigmax] = ...
    fastfindNNPre(A, SubsequenceLength);
[~, ~, ~, ~, meany, ~, sigmay] = ...
    fastfindNNPre(B, SubsequenceLength);

%% compute the matrix profile--P_AB
current_step = floor(SubsequenceLength*step);
pickedIdx = 1:current_step:ASubsequenceCount;

dotproduct = zeros(ASubsequenceCount,1);
% dotproduct = inf(BSubsequenceCount,1);
refine_distance = inf(ASubsequenceCount,1);
% orig_index = 1:ASubsequenceCount;
m = SubsequenceLength;

disp('Running PreSCRIMP-time series join P_AB:');

for i = 1:length(pickedIdx)
    % compute the distance profile
    idx = pickedIdx(i);
    subsequence = B(idx:idx+SubsequenceLength-1);
    [distanceProfile, ~, ~, ~, ~] = fastfindNN(X, subsequence, n, SubsequenceLength, ...
        sumx2, sumx, meanx, sigmax2, sigmax);
    distanceProfile = abs(distanceProfile);
    
    %-------------------added by Sarah----------------------------------
    % replace distance profile correspondings to subsequence contains NaN
    % with inf
    distanceProfile(isSkipA) = inf; 
    %-------------------------------------------------------------------
    
    % figure out and store the nearest neighbor informtion of current step
    % locations (PickedIdx)
    if i == 1
        MatrixProfile = distanceProfile;
        MPindex(:) = idx;
%         [MatrixProfile(idx), MPindex(idx)] = min(distanceProfile);
    else
        updatePos = distanceProfile < MatrixProfile;
        MPindex(updatePos) = idx;
        MatrixProfile(updatePos) = distanceProfile(updatePos);
%         [MatrixProfile(idx), MPindex(idx)] = min(distanceProfile);
    end
    
    idx_nn = MPindex(idx);
%     idx_diff = idx_nn-idx;%--?
    dotproduct(idx) = (m - MatrixProfile(idx)^2/2)*sigmax(idx)*sigmay(idx_nn) + m*meanx(idx)*meany(idx_nn);
    
    t_end_idx_B = min([idx_nn+current_step-1,length(B)]);
    t_end_idx_A = min([idx   +current_step-1,length(A)]);
    fwd_cnt = min([t_end_idx_B-idx_nn+1,t_end_idx_A-idx+1,ASubsequenceCount-idx+1,ASubsequenceCount-idx_nn+1]);

    %endidx = min([ASubsequenceCount, idx+current_step-1, BSubsequenceCount]);%end index of MP(Q). 
%     calculate forward MPs based on current MP location
    % matlab cumsum: B = cumsum(A); If A is a vector, then cumsum(A) returns a vector containing the cumulative sum of the elements of A.
    
    %TODO
%     disp([length(dotproduct),idx+fwd_cnt,length(A),idx+m+fwd_cnt-1,length(B),idx_nn+m+fwd_cnt-1,length(meanx),idx+fwd_cnt,length(meany),idx_nn+fwd_cnt]);
    %TODO
    % 1751        1801        2050        2100        2050        1725        1751        1801        1751        1426

    dotproduct(idx+1:idx+fwd_cnt) =  dotproduct(idx)+...
                                     cumsum( A(idx+m:idx+m+fwd_cnt-1).*B(idx_nn+m:idx_nn+m+fwd_cnt-1) - A(idx:idx+fwd_cnt-1).*B(idx_nn:idx_nn+fwd_cnt-1) );
    refine_distance(idx+1:idx+fwd_cnt) = sqrt(abs(2*(m-(dotproduct(idx+1:idx+fwd_cnt)-m*meanx(idx+1:idx+fwd_cnt).*meany(idx_nn+1:idx_nn+fwd_cnt))./(sigmax(idx+1:idx+fwd_cnt).*sigmay(idx_nn+1:idx_nn+fwd_cnt)))));
    update_pos1 = find(refine_distance(idx+1:idx+fwd_cnt)<MatrixProfile(idx+1:idx+fwd_cnt));
    MatrixProfile(update_pos1) = refine_distance(update_pos1);
    update_j1 = update_pos1 - idx + MPindex(update_pos1);
    MPindex(update_pos1) = update_j1;
    % Need to update MPI
    

    %calculate backward MPs based on current MP location
    t_bg_idx_B = max([idx_nn-current_step+1,1]);
    t_bg_idx_A = max([idx   -current_step+1,1]);
    bwd_cnt = min([idx_nn-t_bg_idx_B+1,idx-t_bg_idx_A+1]);
    %     beginidx = max([1, idx-current_step+1, 1-idx_diff]);%calculate backward MPs based on current MP location
   
    %TODO
%     disp([idx,bwd_cnt,idx-bwd_cnt,idx_nn,idx_nn-bwd_cnt,idx-1+m,idx-bwd_cnt+m,idx_nn-1+m,idx_nn-bwd_cnt+m]);
    %TODO
    if idx > 1 && idx_nn > 1
        dotproduct(idx-1:-1:idx-bwd_cnt) = dotproduct(idx)+...
                                           cumsum(A(idx-1:-1:idx-bwd_cnt).*B(idx_nn-1:-1:idx_nn-bwd_cnt)-A(idx-1+m:-1:idx-bwd_cnt+m).*B(idx_nn-1+m:-1:idx_nn-bwd_cnt+m));
        refine_distance(idx-bwd_cnt:idx-1) = sqrt(abs(2*(m-(dotproduct(idx-bwd_cnt:idx-1)-m*meanx(idx-bwd_cnt:idx-1).*meany(idx_nn-bwd_cnt:idx_nn-1))./(sigmax(idx-bwd_cnt:idx-1).*sigmay(idx_nn-bwd_cnt:idx_nn-1)))));
        update_pos2 = find(refine_distance(idx-1:-1:idx-bwd_cnt)<MatrixProfile(idx-1:-1:idx-bwd_cnt));
        MatrixProfile(update_pos2) = refine_distance(update_pos2);
        
        update_j2 = MPindex(update_pos2) - (idx - update_pos2);
        MPindex(update_pos2) = update_j2;
        % Need to update MPI
        % MPindex(update_pos2+beginidx+idx_diff-1) = orig_index(update_pos2+beginidx+idx_diff-1)-idx_diff;
    end
    
    if mod(i,1000)==0
        completed_perc = i/length(pickedIdx)*100;
        disp(['progress:', num2str(completed_perc), '%']);
    end
end

PreMatrixProfile = MatrixProfile;
PreMPindex = MPindex;

disp('PreSCRIMP-time series join finished');

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