function [G,Av_dist,C_G] = SmallWorld(n,K,p,C_G_min,connected,verbose)
% INPUT
% n is the # of nodes of the generated graph G
% K>0 is an integer s.t. nK is close to the # of edges |Eo| of the und.cted
%   version of the generated graph, if C_G_min = 0! In general: |Eo| >= nK
% p is a probability in [0,1]
% C_G_min is the minimum required glabal clustering coefficient, in [0,1]
% connected is a flag in {0,1}, connected = 1 forces G to be wk. connected
% computes 2nd and 3rd outputs if verbose = 1

% OUTPUT
% G is the generated graph -> Go = (Vo,Eo): undirected version of G
% Av_dist is the average distance in Go (if Go connected: Av_dist ~ log(n))
% C_G >= C_G_min is the glabal clustering coefficient of Go


%% 1 - builds a skeleton through the Watts-Strogatz algorithm
% (the resulting graph has exactly n nodes and nK edges)
if ~(floor(n/2) > K && n > 3) %  && K > log(n)/2
    fprintf('\n')
    fprintf('Please provide n: n > 3 and K: log(n)/2 < K < floor(n/2).\n')
    error('Cannot build a Watts Strogatz skeleton with these data.')
end
A = RegularRingLattice(n,K);    % adjacency matrix of the generated graph
for i = 1:n
    for j = 1:n
        if A(i,j) == 1 && rand <= p
            kk_ij = ceil(rand*(n-sum(A(i,:))-1));
            k = 0;
            while kk_ij > 0
                k = k + 1;
                if A(i,k) == 0 && i ~= k
                    kk_ij = kk_ij - 1;
                    if kk_ij == 0
                        A(i,j) = 0;
                        A(j,i) = 0;
                        A(i,k) = 1;
                        A(k,i) = 1;
                    end
                end
            end
        end
    end
end


%% 2 - makes the skeleton connected if it is not and it is required
if connected
    G = graph(A);
    [bins,binsizes] = conncomp(G);
    NCC = length(binsizes);             % number of connected components
    if NCC > 1
        v_star = zeros(NCC,1);          % nodes to be connected..
        degs_star = zeros(NCC,1)+Inf;   % ..should have lowest degree
        for i = 1:n
            cc = bins(i);               % connected comp. node i belongs to
            if degree(G,i) < degs_star(cc)
                degs_star(cc) = degree(G,i);
                v_star(cc) = i;
            end
        end
        % constructs a random connected subgraph linking all the nodes 
        % in v_star
        for i = 1:NCC
            for j = 1:NCC
                if i ~= j
                    if rand <= (1-p) || i-1 == j
                        A(v_star(i),v_star(j)) = 1;
                        A(v_star(j),v_star(i)) = 1;
                    end
                end
            end
        end
    end
end


%% 3 - assures global clustering coefficient >= C_G_min
C_G = 0;
if C_G_min > 0 && verbose
    C_G = C(A);
end
i = 1;
while C_G < C_G_min && i <= n
    j = 1;
    while C_G < C_G_min && j <= n
        k = 1;
        while C_G < C_G_min && k <= n
            if i ~= j && j ~= k && k ~= i
                cond1 = A(i,j) == 0 && A(j,k) == 1 && A(k,i) == 1;
                cond2 = A(i,j) == 1 && A(j,k) == 0 && A(k,i) == 1;
                cond3 = A(i,j) == 1 && A(j,k) == 1 && A(k,i) == 0;
                if cond1 || cond2 || cond3
                    A(i,j) = 1;
                    A(j,i) = 1;
                    A(j,k) = 1;
                    A(k,j) = 1;
                    A(k,i) = 1;
                    A(i,k) = 1;
                    C_G = C(A);
                end
            end
            k = k + 1;
        end
        j = j + 1;
    end
    i = i + 1;
end


%% 4 - computes average distance
Av_dist = -1;
if verbose
    Av_dist = sum(sum(distances(graph(A))))/(n*(n-1));
end


%% 5 - renders the generated small-world graph directed
for i = 2:n
    for j = 1:i-1
        if A(i,j) == 1
            action = rand;
            if 1/3 <= action && action < 2/3
                A(j,i) = 0;
            elseif 2/3 <= action && action < 1
                A(i,j) = 0;
            end
        end
    end
end


%% 6 - final results (Av_dist and C_G are computed in 4 and 3, resp.)
G = digraph(A);

% figure
% plot(graph(A))
% pause

end


%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% generates a regular ring lattice graph
function A = RegularRingLattice(n,K)

if K < 1 || K >= floor(n/2)
    fprintf('It has to hold that 0 < K < floor(n/2).\n')
    error('Cannot construct a regular lattice graph for this (n,K).\n')
end

w = zeros(1,n);
for i = 1:K
    w(i+1) = 1;
    w(n-i+1) = 1;
end
A = zeros(n);
for i = 1:n
    A(i,:) = circshift(w,[0 i-1]);
end

end


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