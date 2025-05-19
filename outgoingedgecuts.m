function OECs = outgoingedgecuts(A,k,idx)
% INPUT
% A is the adjacency matrix of the given node G
% k is the total number of clusters (it should be k==max(idx))
% idx(i) contains the community to which node i belongs

% OUPUT
% OECs is a cell array containing all the transpose of the adjacency
%    matrices of the outgoing edge-cut subgraphs
%    associated to each cluster according to idx


[~,n] = size(A);
OECs = cell(1,k);
for kk = 1:k

    % characterizing the outgoing edge-cut of cluster kk
    % +1: node inside the cluster kk
    % -1: node outside the cluster kk but reached
    %  0: unreached node outside the cluster kk
    charCL_kk = zeros(1,n);
    for i = 1:n
        if idx(i) == kk
            charCL_kk(i) = 1;
            for j = find(A(i,:))
                if idx(j) ~= kk
                    charCL_kk(j) = -1;
                end
            end
        end
    end

    % finding the interesting nodes for the outgoing edge-cut of cluster kk
    % to be mapped
    n_kk = sum(abs(charCL_kk));
    v_kk = zeros(1,n_kk);
    i_kk = 0;
    for i = 1:n
        if abs(charCL_kk(i))
            i_kk = i_kk + 1;
            v_kk(i_kk) = i;
        end
    end

    % computing the transpose of the adjacency matrix of the 
    % outgoing edge-cut of cluster kk
    OECs{kk} = zeros(n_kk,n_kk);
    for i_kk = 2:n_kk
        i = v_kk(i_kk);
        for j_kk = 1:i_kk-1
            j = v_kk(j_kk);
            if charCL_kk(i) == 1
                OECs{kk}(i_kk,j_kk) = A(i,j);
            end
            if charCL_kk(j) == 1
                OECs{kk}(j_kk,i_kk) = A(j,i);
            end
        end
    end

end

end