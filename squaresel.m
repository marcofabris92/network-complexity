%% selects rows and columns of matrix M according to components in v
function M_v = squaresel(M,v)

l = length(v);
M_v = zeros(l,l);
for i = 1:l
    i_v = v(i);
    M_v(i,i) = M(i_v,i_v);
    for j = 1:i-1
        j_v = v(j);
        M_v(i,j) = M(i_v,j_v);
        M_v(j,i) = M(j_v,i_v);
    end
end

end