function G = ErdosRenyi(n,p)
% INPUT
% n is the # of nodes of the generated graph
% p is a probability in [0,1]

% OUTPUT
% G is the generated graph


%% 1 - generation of the graph through Erdos-Renyi algorithm
G = zeros(n,n);     % adjacency matrix of the generated graph
for i = 1:n
    for j = 1:n
        if i ~= j && rand <= p
            G(i,j) = 1;
        end
    end
end


%% 2 - final result
G = digraph(G); 

end