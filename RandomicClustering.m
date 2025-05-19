function [idx,idxsizes] = RandomicClustering(G,k,type,verbose)
% INPUT
% G is the given graph
% k is the MAXIMUM number of clusters that this algorithm can make
% type = 'random'
%        'inrandom'
%        'outrandom'
%        'neutralrandom'
%        'polarizedrandom'
%        'purelyrandom'
%        'naive'
% verbose = 1 -> computes 2nd output

% OUPUT
% idx(i) contains the community to which node i belongs
% idxsizes(j) tells the size of cluster j


%% 0A - naive clustering
if strcmp(type,'naive')
    n = numnodes(G);
    idx = zeros(1,n);
    idxsizes = zeros(1,k);
    ratio = floor(n/k);
    remainder = mod(n,k);
    i = 0;
    while i < k
        i = i + 1;
        for j = 1:ratio
            idx((i-1)*ratio+j) = i;
        end
        if verbose
            idxsizes(i) = ratio;
        end
        if i <= remainder
            idx(k*ratio+i) = i;
            if verbose
                idxsizes(i) = idxsizes(i) + 1;
            end
        end
    end
end


%% 0B - purerandom clustering
if strcmp(type,'purelyrandom')
    n = numnodes(G);
    idx = ones(1,n);
    idxsizes = k;
    if k == n
        idx = 1:n;
        if verbose
            idxsizes = ones(1,k);
        end
    elseif 1 < k && k < n
        idx = ceil(k*rand(1,n));
        idxsizes = clcard(idx,k);
        for kk = 1:k
            if idxsizes(kk) == 0
                [maxsize,kk_maxsize] = max(idxsizes);
                direction = floor(2*rand);
                i = 1;
                if direction
                    i = n;
                end
                des_size = floor(maxsize/2);
                while idxsizes(kk) < des_size
                    if idx(i) == kk_maxsize
                        idx(i) = kk;
                        idxsizes(kk) = idxsizes(kk) + 1;
                        idxsizes(kk_maxsize) = idxsizes(kk_maxsize) - 1;
                    end
                    if direction
                        i = i - 1;
                    else
                        i = i + 1;
                    end
                end
            end
        end
    end
    return
end



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


%% 2C_a - random clustering
c_scl = 0;    % counter of subclusters
p = 1;
for cc = 1:NCC
    K_cc = K(cc);
    [nodes_cc,n_cc] = nodesinacc(bins,binsizes,cc);
    
    idx_cc = ones(n_cc,1);
    if K_cc > 1
        dIN = 0;
        dOUT = 0;
        dd = 0;

        if strcmp(type,'inrandom') || strcmp(type,'outrandom') ||... 
           strcmp(type,'neutralrandom') || strcmp(type,'symmetrizedrandom')
            for cc2 = 1:n_cc
                dcc2_in = indegree(G,nodes_cc(cc2));
                dcc2_out = outdegree(G,nodes_cc(cc2));
                dIN = dIN + dcc2_in;
                dOUT = dOUT + dcc2_out;
                dd = dd + max(dcc2_in,dcc2_out);
            end
        end

        for cc2 = 1:n_cc

            if strcmp(type,'inrandom')
                dcc2 = indegree(G,nodes_cc(cc2));
                p = (1+dcc2)/(1+dIN);
            end
            if strcmp(type,'outrandom')
                dcc2 = outdegree(G,nodes_cc(cc2));
                p = (1+dcc2)/(1+dOUT);
            end
            if strcmp(type,'neutralrandom')
                dcc2_in = indegree(G,nodes_cc(cc2));
                dcc2_out = outdegree(G,nodes_cc(cc2));
                dcc2 = max(dcc2_in,dcc2_out);
                p = (dd/(dIN+dOUT)+dcc2)/(dd/(dIN+dOUT)+dd);
            end
            if strcmp(type,'symmetrizedrandom')
                dcc2_in = indegree(G,nodes_cc(cc2));
                dcc2_out = outdegree(G,nodes_cc(cc2));
                p = (2+dcc2_in+dcc2_out)/(2+dIN+dOUT);
            end

            if cc2 <= K_cc
                idx_cc(cc2) = cc2;
            else
                if rand < p
                    idx_cc(cc2) = ceil(rand*K_cc);
                else
                    idx_cc(cc2) = ceil(p*K_cc);
                end
            end

        end
    end

    for i = 1:n_cc
        idx(nodes_cc(i)) = idx_cc(i)+c_scl;
    end
    c_scl = c_scl + K_cc;
end


%% 2C - completion
if verbose
    % check_ks = zeros(1,k);
    % for kk = 1:k
    %     check_ks(idx(kk)) = 1;
    % end
    % 
    % jj = Inf;
    % true_k = k;
    % while jj > 0
    %     kk = 0;
    %     jj = 0;
    %     while kk < true_k && jj == 0
    %         kk = kk + 1;
    %         if check_ks(kk) == 0
    %             check_ks(kk) = 1;
    %             jj = kk;
    %         end
    %     end
    %     if jj > 0
    %         true_k = true_k - 1;
    %         [maxidx,posj] = max(idx);
    %         check_ks(maxidx) = 0;
    %         for j = 1:length(posj)
    %             idx(posj(j)) = jj;
    %         end
    %     end
    % end
    % idx = idx - min(idx) + 1;
    idxsizes = clcard(idx,k);
end


end

