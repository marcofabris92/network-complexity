function [bins,binsizes,NCC,idx,idxsizes] = initclustering(G,verbose)

[bins,binsizes] = conncomp(G,'Type','weak');
NCC = length(binsizes);
idx = bins;
idxsizes = 0;
if verbose
    idxsizes = clcard(idx,NCC);
end

end