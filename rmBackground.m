function new_trees = rmBackground(trees, fmap, params)
% RMBACKGROUND  Remove every first level subtrees which have no node located on foreground

[rows cols nfms] = size(fmap);
N = rows * cols;

for i = 1:length(trees)
    tree = trees{i};
    P = [tree(:).Parent];
    idx = find(P(2:end) == 1)+1; % the root is the whole frame
    dt = treeDecendt(tree);
    nt = length(tree);
    S = zeros(nt, 1);
    offset = (i - 1) * N;
    for j = 1:nt
        S(j) = sum(fmap(offset + tree(j).PixelIdxList))/ length(tree(j).PixelIdxList);
    end
    tree1 = tree(1);
    for j = 1:length(idx)
        id = idx(j);
        idx2 = [id dt{id}];
        if any(S(idx2) > params.er_score_thre)
            tree1 = [tree1 tree(idx2)];
        end    
    end
    new_trees{i} = tree1;
end