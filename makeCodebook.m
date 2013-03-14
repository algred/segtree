function [idx C] = makeCodebook(X, prm)

[L, C, d] = kmeans2(X, prm.k, prm );

c = unique(L);
D = pdist2(X, C, prm.metric);
for i = 1:length(c)
    idx1 = find(L == c(i));
    [~, id] = min(D(idx1, i));
    idx(i) = idx1(id);
end
C = X(idx, :);

end