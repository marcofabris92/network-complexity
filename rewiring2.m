% rewiring preserving number of nodes, number of edges and weights
% tends to create nodes with high outdegree if p -> 1
% boils down to rewiring1 if p -> 0
function G = rewiring2(G,par)

A = zeros(par.n,par.n);
k = 0;
while true
    rng("shuffle")
    I = 1:par.n;
    J = I(randperm(par.n));
    I = I(randperm(par.n));
    for i = I
        for j = J
            if i ~= j && A(i,j) == 0 && rand < par.p
                k = k + 1;
                A(i,j) = par.weights(k);
                if k == par.m
                    G = digraph(A);
                    %[numnodes(G) numedges(G)]
                    return
                end
            end
        end
    end
end

end