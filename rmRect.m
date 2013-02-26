function tree = rmRect(tree, rows, cols, params)

mask = zeros(rows, cols);
keep = ones(size(tree));
P = [tree(:).Parent];
sidx = find(P(2:end) == 1) + 1;
D = treeDecendt(tree);
d = params.er_delta;
for i = 1:length(sidx)
    id = sidx(i);
    mask(tree(id).PixelIdxList) = 1;
    B = bwboundaries(mask);
    x = B{1}(:,2);
    y = B{1}(:,1);
    dx1 = x - x([d:end 1:(d-1)]);
    dx2 = x([end-d+1:end 1:end-d]) - x;
    dy1 = y - y([d:end 1:(d-1)]);
    dy2 = y([end-d+1:end 1:end-d]) - y;
    vflg = (dy1 == 0 & dy2 == 0);
    hflg = (dx1 == 0 & dx2 == 0);
    flg = ~(vflg | hflg);
    lflg = (abs(dx1(flg) ./ (dy1(flg)+eps) - dx2(flg) ./ (dy2(flg) + eps)) < 0.001);
    r = (sum(vflg | hflg) + sum(lflg)) / (length(x) - 2);
    if r > 0.5
        keep([id D{id}]) = 0;
    end
    mask(tree(id).PixelIdxList) = 0;
end

idx = find(keep);
for i = 2:length(idx)
    p2 = find(P == idx(i));
    for j = 1:length(p2)
        tree(p2(j)).Parent = i;
    end
end
tree = tree(keep > 0);

end