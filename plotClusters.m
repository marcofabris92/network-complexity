function [] = plotClusters(G,k,idx,idxsizes)
% INPUT
% G is the given graph
% k is the total number of clusters
% idx(i) tells which cluster node i belongs to
% idxsizes(j) tells the size of cluster j

%% OUPUT
% depicts clusters in G according to idx

figure('Position',[200 200 600 600])
h = plot(G);
hold on
h.EdgeColor = [0 0 0];
h.EdgeAlpha = 1;
h.MarkerSize = 7;
color = hsv(k);

for kk = 1:k
    hl_kk = zeros(1,idxsizes(kk));
    j = 0;
    for i = 1:numnodes(G)
        if idx(i) == kk
            j = j + 1;
            hl_kk(j) = i;
        end
    end
    highlight(h,hl_kk,'NodeColor',color(kk,:))
end

end
