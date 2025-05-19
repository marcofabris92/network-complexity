function [idx,idxsizes] = underclustering(k,bins,binsizes,idx,verbose)

CL = kmeans(binsizes',k);
for i = 1:length(bins)
    idx(i) = CL(bins(i));
end
idxsizes = 0;
if verbose
    idxsizes = clcard(idx,k);
end

end
