function tree = treeDepth(tree)

P = [tree(:).Parent];
depth = ones(size(P));

for i = 2:length(P)
    depth(i) = depth(P(i)) + 1;
end

for i = 1:length(tree)
    tree(i).Depth = depth(i);
end

end