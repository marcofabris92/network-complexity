function G = rewiring(G,par,type)

switch type 
    case 0
        G = rewiring0(G,par);
    case 1
        G = rewiring1(G,par);
    case 2
        G = rewiring2(G,par);
    case 3
        G = rewiring3(G,par);
    otherwise
        error('incorrect type')
end

end