function [idx,idxsizes] = CentralityClustering(G,k,type,verbose)
% INPUT
% G is the given graph
% k is the fix desired number of clusters
% type = 'indegree'
%        'outdegree'
%        'incloseness'
%        'outcloseness'
%        'betweenness'
%        'pagerank'
%        'hubs'
%        'authorities'

% OUTPUT
% idx(i) contains the community to which node i belongs
% idxsizes(j) tells the size of cluster j


kmeans_warn = 'stats:kmeans:FailedToConverge';
warning('off',kmeans_warn);
idx = kmeans(centrality(G,type),k);
warning('on',kmeans_warn);
if verbose
    idxsizes = clcard(idx,k);
end

    
end