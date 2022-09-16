function Z = pick_idx(cluster_size, N)
% given a vector of cluster sizes, return the corresponding indices the
% vectorized edge list.

temp = zeros(N);

temp(1:cluster_size,1:cluster_size) = 1;

for j = 1:N
    temp(j,j) = 0;
end
Z= find(squareform(temp));

end