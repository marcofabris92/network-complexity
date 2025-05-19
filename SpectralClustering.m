function [idx,idxsizes] = SpectralClustering(G,k,type,verbose)
% INPUT
% G is the given graph
% k is the MAXIMUM number of clusters that this algorithm can make
% type = 'neutral'
%        'neutralsquared'
%        'symmetrized'
%        'symmetrizedsquared'
%        'inward'
%        'outward'
%        'inoutward'
%        'outinward'
%        'strong'

% OUPUT
% idx(i) contains the community to which node i belongs
% idxsizes(j) tells the size of cluster j


%% 1 - computes weakly connected components
[bins,binsizes,NCC,idx,idxsizes] = initclustering(G,verbose);


%% 2A - Component clustering, if k == NCC
if k == NCC
    return
end


%% 2B - Component clustering, if k < NCC
if k < NCC
    [idx,idxsizes] = underclustering(k,bins,binsizes,idx,verbose);
    return
end


%% 2C - Component clustering, if k > NCC
[K,bins,binsizes] = subcldims(k,bins,binsizes);


%% 2C_a - various types of spectral clustering (not strong)
if strcmp(type,'neutral') || strcmp(type,'symmetrized') || ...
      strcmp(type,'inward') || strcmp(type,'outward') ||...
      strcmp(type,'neutralsquared') || strcmp(type,'symmetrizedsquared')...
      || strcmp(type,'inoutward') || strcmp(type,'outinward')

    L = full(adjacency(G,'weighted'));
    if strcmp(type,'neutral')
        L = undLap(L);
    elseif strcmp(type,'neutralsquared')
        L = undLap(L);
        L = L*L;
    elseif strcmp(type,'symmetrized')
        L = diag(sum(L+L'))-(L+L');
    elseif strcmp(type,'symmetrizedsquared')
        L = diag(sum(L+L'))-(L+L');
        L = L'*L;
    elseif strcmp(type,'inward')
        L = diag(sum(L,1))-L;
        L = L'*L;
    elseif strcmp(type,'outward')
        L = diag(sum(L,2))-L;
        L = L'*L;
    elseif strcmp(type,'inoutward')
        LI = diag(sum(L,1))-L;
        LO = diag(sum(L,2))-L;
        L = LI'*(LO'*LO)*LI;
    elseif strcmp(type,'outinward')
        LI = diag(sum(L,1))-L;
        LO = diag(sum(L,2))-L;
        L = LO'*(LI'*LI)*LO;
    end

    c_scl = 0;    % counter of subclusters
    for cc = 1:NCC
        K_cc = K(cc);
        [nodes_cc,n_cc] = nodesinacc(bins,binsizes,cc);
        
        idx_cc = ones(n_cc,1);
        if K_cc > 1
            s = -1;
            if strcmp(type,'neutral') || strcmp(type,'neutralsquared')...
                    || strcmp(type,'symmetrized') ||...
                    strcmp(type,'symmetrizedsquared')
                s = 1;
            elseif strcmp(type,'inward') || strcmp(type,'outward') ||...
                    strcmp(type,'inoutward') || strcmp(type,'outinward') 
                s = 0;
            end
            idx_cc = spectralkmeans(K_cc,squaresel(L,nodes_cc),s);
        end

        for i = 1:n_cc
            idx(nodes_cc(i)) = idx_cc(i)+c_scl;
        end
        c_scl = c_scl + K_cc;
    end
end


%% 2C_b - clustering based on strongly connected components
if strcmp(type,'strong')

    A = adjacency(G);

    c_scl = 0;
    for cc = 1:NCC
        K_cc = K(cc);
        [nodes_cc,n_cc] = nodesinacc(bins,binsizes,cc);

        Acc = sparse(squaresel(A,nodes_cc));
        [bins_cc,binsizes_cc] = conncomp(digraph(Acc),'Type','strong');
        max_bs = max(binsizes_cc);
        for i = 1:n_cc
            b_cc_i = bins_cc(i);
            bs_cc_i = binsizes_cc(b_cc_i);
            for j = 1:n_cc
                b_cc_j = bins_cc(j);
                bs_cc_j = binsizes_cc(b_cc_j);
                if i ~= j && b_cc_i ~= b_cc_j && Acc(i,j) > 0
                    Acc(i,j) = ((bs_cc_i+bs_cc_j)/2-1)/max_bs;
                end
            end
        end
        idx_cc = SpectralClustering(digraph(Acc),K_cc,...
            'symmetrized',verbose);

        for i = 1:n_cc
            idx(nodes_cc(i)) = idx_cc(i)+c_scl;
        end
        c_scl = c_scl + K_cc;
    end
end


%% 2C - completion
if verbose
    idxsizes = clcard(idx,k);
end


end




%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% spectral clustering + kmeans executed on the Laplacian L_cc
function idx_cc = spectralkmeans(K_cc,L_cc,s)

[V,D] = eig(L_cc);
[~,ind] = sort(diag(D),'ComparisonMethod','abs');
V = V(:,ind);
idx_cc = kmeans(V(:,1+s:K_cc),K_cc,'start',zeros(K_cc,K_cc-s));

end


%% retrieves undirected Laplacian from a directed adjacency matrix
function L = undLap(A)

n = length(A(1,:));
L = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j && A(i,j) > 0
            if A(j,i) == 0
                L(i,j) = A(i,j);
                L(j,i) = L(i,j);
            else
                av = mean([A(i,j) A(j,i)]);
                L(i,j) = av;
                L(j,i) = av;
            end
        end
    end
end
L = diag(sum(L))-L;

end

