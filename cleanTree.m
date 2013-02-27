function tree = cleanTree(tree, params)
% CLEANTREE    Remove redundant segments from a segment tree


P = [tree(:).Parent];
keep = ones(1, length(tree));

N = length(tree(1).PixelIdxList);
maxsz = params.maxsz_ratio * N;

for i = 2:length(tree)
    sz = length(tree(i).PixelIdxList);
    sz1 = length(tree( tree(i).Parent ).PixelIdxList);
    if sz >= sz1 * params.min_shrink_ratio || sz >= maxsz
        child = find( P == i );
        for j = 1:length(child)
            tree( child(j) ).Parent = tree(i).Parent;
        end
        P( P == i ) = tree(i).Parent;
        keep(i) = 0;
    end
end

idx = find(keep > 0);
idmap = containers.Map(idx, 1:length(idx));
tree = tree(idx);
for i = 1:length(tree)
    tree(i).Parent = idmap(tree(i).Parent);
end

end