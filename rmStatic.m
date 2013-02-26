function tree = rmStatic(tree, M, params)

[rows cols] = size(M);
N = rows * cols;

% estimate background motion threshold
bgidx = [];
for i = 2:length(tree)
    tree(i).Old_ID = i;
    bgidx = union(bgidx, tree(i).PixelIdxList);
end
bgidx = setdiff(1:N, bgidx);
bgm = round(M(bgidx)*10);
maxM = max(bgm);
bgp = cumsum(histc(bgm, 0:maxM))/length(bgidx);
tmp = find(bgp > 0.5);
thre = tmp(1)/10;

% remove first level subtrees in which no node have average motion
% magnitude larger than the threshold
m = zeros(size(tree));
for i = 2:length(tree)
    m(i) = sum(M(tree(i).PixelIdxList)) / length(tree(i).PixelIdxList);
end
m = m > thre;

P = [tree(:).Parent];
stidx = find(P(2:end) == 1) + 1;
descendt = treeDecendt(tree);
keep = ones(1, length(tree));
for i = 1:length(stidx)
    if ~any(m(descendt{stidx(i)}))
        keep(descendt{stidx(i)}) = 0;
    end
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