function Phi_G = netcompl_(G,K,clusters)
% netcompl function variant: for manually selected clusters
% see also netcompl.m

Phi_G = 0;
OECs = outgoingedgecuts(adjacency(G),K,clusters);
for k = 1:K
    Phi_G = Phi_G + BipGraph(OECs{k}).hopcroftKarp();
end

end