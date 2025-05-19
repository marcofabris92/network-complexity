% rewiring preserving number of nodes, number of edges and weights
% tends to create balanced digraphs
function G = rewiring1(G,par)

A = zeros(par.n,par.n);
k = 0;
while true
    rng("shuffle")
    for q = par.Q(randperm(length(par.Q)))
        j = ceil(q/(par.n-1));
        i = q-(j-1)*(par.n-1);
        i = i + (i >= j);
        if A(i,j) == 0 && rand < par.p
            k = k + 1;
            A(i,j) = par.weights(k);
            if k == par.m
                G = digraph(A);
                return
            end
        end
    end
end

end