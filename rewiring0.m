% rewiring preserving number of nodes, number of edges, weights and outdegrees
function G = rewiring0(G,par)

A = adjacency(G);
for i = 1:par.n
    for j = [1:i-1 i+1:par.n]
        if A(i,j) ~= 0 && rand < par.p
            k = randi(par.n-1);
            k = k + (k >= i);
            G = addedge(rmedge(G,i,j),i,k,A(i,j));
        end
    end
end

end