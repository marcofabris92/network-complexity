% BipGraph class for graph transformation and Hopcroft-Karp algorithm
classdef BipGraph < handle

    properties (SetAccess = 'private', GetAccess = 'private')
        n
        Adj
        Pair_U
        Pair_V
        Dist
        max_Odeg
    end

    %% methods
    % BipGraph          % constructor
    % hopcroftKarp      % returns size of maximum matching
    % bfs               % breadth-first search (= true if augmenting path
                        % it is found)
    % dfs               % depth-first search (adds augmenting path if there
                        % is one beginning with u
    methods

        %% constructor
        function obj = BipGraph(A)
            obj.n = length(A(1,:));
            obj.Pair_U = NaN*ones(obj.n,1);
            obj.Pair_V = NaN*ones(obj.n,1);
            obj.Dist = -ones(obj.n+1,1);
            obj.Adj = cell(obj.n,1);
            obj.max_Odeg = 0;
            for i = 1:obj.n
                Odeg_i = sum(A(i,:));
                if Odeg_i > obj.max_Odeg
                    obj.max_Odeg = Odeg_i;
                end
                obj.Adj{i} = zeros(1,Odeg_i);
                j_ = 0;
                for j = 1:obj.n
                    if A(i,j) == 1
                        j_ = j_ + 1;
                        obj.Adj{i}(j_) = j;
                    end
                end
            end
        end

        %% hopcroftKarp
        function matching_size = hopcroftKarp(obj)
            matching_size = 0;
            while obj.bfs()
                for u = 1:obj.n
                    if isnan(obj.Pair_U(u)) && obj.dfs(u)
                        matching_size = matching_size + 1;
                    end
                end
            end
        end

    end

    methods (Access = private)

        %% bfs
        function found = bfs(obj)
            Q = FiniteQueue(obj.max_Odeg+1);
            for u = 1:obj.n
                if isnan(obj.Pair_U(u))
                    obj.Dist(u) = 0;
                    Q.enqueue(u);
                else
                    obj.Dist(u) = +Inf;
                end
            end
            obj.Dist(end) = +Inf;
            while ~Q.isEmpty()
                u = Q.front();
                Q.dequeue();
                if ~isnan(u) && obj.Dist(u) < obj.Dist(end)
                    for v = obj.Adj{u}
                        u_prime = obj.Pair_V(v);
                        if isnan(u_prime) && ...
                               obj.Dist(end) == +Inf || ...
                               ~isnan(u_prime) && ...
                               obj.Dist(u_prime) == +Inf
                            if isnan(u_prime)
                                obj.Dist(end) = obj.Dist(u)+1;
                            else
                                obj.Dist(u_prime) = obj.Dist(u)+1;
                            end
                            Q.enqueue(u_prime);
                        end
                    end
                end
            end
            found = (obj.Dist(end) ~= +Inf);
            Q.del();
        end
        
        %% dfs
        function found = dfs(obj,u)
            if ~isnan(u)
                for v = obj.Adj{u}
                    u_prime = obj.Pair_V(v);
                    if isnan(u_prime) && ...
                            obj.Dist(end) == obj.Dist(u)+1 || ...
                            ~isnan(u_prime) && ...
                            obj.Dist(u_prime) == obj.Dist(u)+1
                        if obj.dfs(u_prime)
                            obj.Pair_V(v) = u;
                            obj.Pair_U(u) = v;
                            found = 1;
                            return
                        end
                    end
                end
                obj.Dist(u) = +Inf;
                found = 0;
                return
            end
            found = 1;
        end

    end
end