function [idx,idxsizes,best_type,best_modularity] =...
    clustering(G,k,verbose,type,dic,dic_info)
% INPUT
% G is the input graph
% k is number of clusters to make
% verbose = 1 -> plots some clusters and computes idxsizes
% type is the clustering method to be used
% dic is a dictionary containing all clustering types
% dic_info specifies extra data for dic

% OUTPUT
% idx(i) contains the community to which node i belongs
% idxsizes(j) tells the size of cluster j
% best_type is the clustering method optimizing modularity
% best_modularity is the optimized modularity value


%% 1 - Checks k and the given type; retrieves dic information
[key_type,k] = checks(G,k,type,dic,dic_info(2:3));
mymod = dic_info(1);
shifted_dim = dic_info(3);
key_trans = dic_info(4:end);


%% 2 - applies the preferred clusterization
idxsizes = 0;
best_type = NaN;
best_modularity = NaN;
if key_type <= 0
    idx_type = zeros(shifted_dim,numnodes(G));
    for key = 1:shifted_dim
        idx_type(key,:) = clustering_(G,k,dic(key),key,key_trans,verbose);
    end
    best_key = 0;
    err_mod = +Inf;
    best_modularity = -Inf;
    if key_type == -1
        best_modularity = +Inf;
    end
    A = adjacency(G);
    for key = 1:shifted_dim
        current_modularity = modularity(A,idx_type(key,:));
        if key_type == 0 && current_modularity > best_modularity ||...
                key_type == -1 && current_modularity < best_modularity
            best_modularity = current_modularity;
            best_key = key;
        end
        if key_type == -2 && abs(current_modularity-mymod) < err_mod
            err_mod = abs(current_modularity-mymod);
            best_modularity = current_modularity;
            best_key = key;
        end
    end
    if best_key == 0 % if modularity has no meaning...
        best_key = getKey(dic,'symmetrized');
    end
    best_type = dic(best_key);
    idx = idx_type(best_key,:);
    if verbose
        idxsizes = clcard(idx,k);
        plotClusters_(G,k,idx,idxsizes);
    end
    return
end

idx = clustering_(G,k,type,key_type,key_trans,verbose);
if verbose
    idxsizes = clcard(idx,k);
    plotClusters_(G,k,verbose,idx,idxsizes)
end


end




%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Checks k and the given type
function [key_type,k] = checks(G,k,type,dic,shifting)

shift = shifting(1);
shifted_dim = shifting(2);
if k < 1 && k > numnodes(G)
    error('It must be 1 <= k <= n.')
end
k = round(k);
key_type = getKey(dic,type);
if isempty(key_type)
    fprintf('The possible options are:\n')
    for key = 1-shift:shifted_dim
        fprintf(['''' dic(key) '''' '\n'])
    end
    error('An invalid option was selected.')
end

end


%% plots clusters
function [] = plotClusters_(G,k,idx,idxsizes)

n = numnodes(G);
if 1 < k && k < n &&...
        (n <= 10 ||...
         n <= 100 && mod(k,10) == 0 ||...
         n <= 1000 && mod(k,100) == 0)
    plotClusters(G,k,idx,idxsizes)
end

end


%% stub function for actual clustering execution
function idx = clustering_(G,k,type,key_type,key_trans,verbose)

if key_type <= key_trans(1)
    idx = CentralityClustering(G,k,type,verbose);
elseif key_type <= key_trans(2)
    idx = SpectralClustering(G,k,type,verbose);
elseif key_type <= key_trans(3)
    idx = RandomicClustering(G,k,type,verbose);
end

end


%% obtains a key from a given value of a dictionary
function key = getKey(dic,value)

testind = cellfun(@(x)isequal(x,value),values(dic));
testkeys = keys(dic);
key = cell2mat(testkeys(testind));

end

