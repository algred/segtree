function cluster = matrix2cluster(M)
% using union find algorithm
n = size(M, 1);
p = 1:n;
idx = find(M(:));
[ys xs] = ind2sub([rows cols], idx);
for i = 1:length(ys)
    p = find(p, ys(i));
    p(xs(i)) = p(ys(i));
end

end


