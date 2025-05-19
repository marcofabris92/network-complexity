%% finds nodes in component cc along with their cardinality n_cc
function [nodes_cc,n_cc] = nodesinacc(bins,binsizes,cc)

n_cc = binsizes(cc);
nodes_cc = zeros(1,n_cc);
j = 0;
for i = 1:length(bins)
    if bins(i) == cc
        j = j + 1;
        nodes_cc(j) = i;
    end
end

end