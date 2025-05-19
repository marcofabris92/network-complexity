% rewiring preserving number of nodes, number of edges and weights
% similarly to rewiring1 but this version is focused on cluster interedges 
function G = rewiring3(G,par)

[sorted_cl,idx] = sort(par.clusters,'ascend');
transitions = zeros(1,par.K);
P = zeros(par.n,par.n);
k = 1;
transitions(k) = 1;
for i = 1:par.n
    if sorted_cl(i) > k
        k = k + 1;
        transitions(k) = i;
    end
    P(i,idx(i)) = 1;
end
A_ = P*adjacency(G)*P';

for k = 1:par.K
    t = transitions(k):(par.n);
    if k < par.K
        t = transitions(k):(transitions(k+1)-1);
    end
    A_k = A_(t,t);
    G_k = digraph(A_k);
    par_k.n = numnodes(G_k);
    par_k.m = numedges(G_k);
    par_k.weights = reshape(A_k,par_k.n^2,1);
    par_k.weights = par_k.weights(par_k.weights~=0);
    par_k.p = par.p;
    par_k.Q = 1:(par_k.n^2 - par_k.n);
    A_(t,t) = adjacency(rewiring1(G_k,par_k));
end

G = digraph(P'*A_*P);

end