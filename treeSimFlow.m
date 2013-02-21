function S = treeSimFlow(tree_prev, tree_cur, u, v)

[rows cols] = size(u);
u1 = floor(u);
v1 = floor(v);
u2 = ceil(u);
v2 = ceil(v);

for i = 1:length(tree_prev)
    sz1(i) = length(tree_prev(i).PixelIdxList);
end

for i = 1:length(tree_cur)
    sz2(i) = length(tree_cur(i).PixelIdxList);
end

for i = 1:length(tree_prev)
    p = tree_prev(i).PixelIdxList;
    xoff1 = u1(p);
    xoff2 = u2(p);
    yoff1 = cols * v1(p);
    yoff2 = cols * v2(p);
    p1 = p + xoff1 + yoff1;
    p2 = p + xoff1 + yoff2;
    p3 = p + xoff2 + yoff1;
    p4 = p + xoff2 + yoff2;
    p = unique([p1 p2 p3 p4]);
    
    for j = 1:length(tree_cur)
        s = length(intersect(p, tree_cur(j).PixelIdxList));
        S(i, j) = s / max(sz1(i), sz2(j));
    end
end