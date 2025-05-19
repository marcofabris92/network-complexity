function idxsizes = clcard(idx,k)
% INPUT
% idx(i) tells which cluster node i belongs to
% k is the total number of clusters

% OUTPUT
% idxsizes(j) tells the size of cluster j


%% counts up the cardinality of each formed cluster
idxsizes = zeros(1,k);
for i = 1:length(idx)
    j = idx(i);
    idxsizes(j) = idxsizes(j) + 1;
end

end