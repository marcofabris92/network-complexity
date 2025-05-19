function Phi_G = netcompl(G,k,verbose,type,dic,dic_info)
% INPUT
% G is the input graph
% k is number of clusters to make
% verbose = 1 -> plots some clusters
% type is the clustering method to be used
% dic is a dictionary containing all clustering types
% dic_info specifies extra data for dic

% OUTPUT
% Phi_G is the network complexity of G


%% Computing Phi_G
Phi_G = 0;
OECs = outgoingedgecuts(adjacency(G),k,...
    clustering(G,k,verbose,type,dic,dic_info));
for kk = 1:k
    Phi_G = Phi_G + BipGraph(OECs{kk}).hopcroftKarp();
end


end