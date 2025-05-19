%% determines subclusters dimensions K and makes order among connected
% components (sorted according to descend size)
function [K,bins,binsizes] = subcldims(k,bins,binsizes)

% Sorts weakly connected components
n = length(bins);
[binsizes,binsidx] = sort(binsizes,'descend');
bi_prev = 0;
j = 0;
for i = 1:n
    bi = bins(i);
    if bi_prev ~= bi
        j = 1;
    end
    bi_prev = bi;
    while binsidx(j) ~= bi 
        j = j + 1;
    end
    bins(i) = j;
end

% creates an array K of subcluster dimensions
NCC = length(binsizes);
K = zeros(1,NCC);
for i = 1:NCC
    K(i) = ceil(k*binsizes(i)/n);
end

% rounds the subcluster dimensions and make them uniform
j = 1;
jlim = NCC+1;
skip = 0;
while skip || sum(K) > k && jlim > 1
    if j == jlim-1 || K(j) > K(j+1)
        skip = 0;
        K(j) = K(j) - 1;
        if K(j) == 1
            jlim = j;
        end
        if j > 1
            j = j - 1;
        end
    else
        skip = 1;
        j = j + 1;
    end
end
if sum(K) ~= k
    error('Something wrong with subcldims function.')
end

end