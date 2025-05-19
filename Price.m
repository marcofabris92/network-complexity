function [G,Av_dist,C_G] = Price(n_ini,m_ini,n,dnew,connected,verbose)
% INPUT
% n_ini is the # of nodes of the initial graph G_ini, s.t. n_ini >= 1
% m_ini >= 0 is the # of edges of the initial graph G_ini
%   if m_ini < n_ini && connected == 1 then m_ini is ignored
% dnew is the out-degree selected for all new nodes, s.t. dnew<=n_ini
% n is the # of nodes of the generated graph G, s.t. n>=n_ini
% connected is a flag in {0,1}, connected = 1 forces G to be wk. connected
% computes 2nd and 3rd outputs if verbose = 1

% OUTPUT
% G is the generated graph -> Go = (Vo,Eo): undirected version of G
% Av_dist is the average distance in Go (if Go connected: Av_dist ~ log(n))
% C_G >= C_G_min is the glabal clustering coefficient of Go


%% 1 - builds the initial graph G_ini
if m_ini < 0 || m_ini > n_ini*(n_ini-1)
    error('It is required 0 <= m_ini <= n_ini*(n_ini-1).')
end
if n_ini < 1 || n_ini > n
    error('It is required 1 <= n_ini <= n.')
end
if dnew > n_ini
    error('It is required dnew <= n_ini.')
end
nt = n_ini;             % dynamic index: tracks the growing number of nodes
A = zeros(n,n);         % adjacency matrix of the generated graph
m_0 = 0;                % counter for initial number of edges
if connected
    for i = 1:n_ini-1
        A(i+1,i) = 1;
    end
    A(1,n_ini) = 1;
    m_0 = n_ini;
end
while m_0 < m_ini
    i = ceil(rand*nt);
    j = ceil(rand*nt);
    if i ~= j && A(i,j) == 0
        A(i,j) = 1;
        m_0 = m_0 + 1;
    end
end


%% 2 - generation of the graph through Price algorithm
% ROW = OUT
% COLUMN = IN
while nt < n
    i = nt+1;                   % new node
    INvol = sum(sum(A));        % current IN-volume of the graph
    C_i = (1:nt);               % canidates out-neighbors set of i
    s_i = 0;                    % counter of selected out-neighbors of i
    while s_i < dnew
        j = C_i(ceil(rand*length(C_i)));    % candidate node to be attached
        if rand <= (1+sum(A(:,j)))/(nt+INvol)
            A(i,j) = 1;
            C_i = setdiff(C_i,j);
            s_i = s_i + 1;
        end
    end
    nt = nt + 1;
end


%% 3 - final results
G = digraph(A);  
Av_dist = -1;
C_G = -1;
if verbose
    Ao = (A+A')/2;
    Av_dist = sum(sum(distances(graph(Ao))))/(n*(n-1));
    C_G = C(Ao);
end

end


%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% global clustering coefficient
function value = C(A)

degs = sum(A);
n = length(degs);
denC = sum(degs.*(degs-1));
value = 0;
for i = 1:n
    for j = 1:n
        for k = 1:n
            value = value + A(i,j)*A(j,k)*A(k,i);
        end
    end
end
value = value/denC;

end