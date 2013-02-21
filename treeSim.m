function S = treeSim(tree_prev, tree_cur, chs, params)
% for a pair of segments, compute its similarity by appearance features

S = zeros(length(tree_prev), length(tree_cur));

for i = 2:length(tree_prev)
    for j = 2:length(tree_cur)
        S(i, j) = getSimilarity(tree_prev(i), tree_cur(j), chs, params.cov_disturb, params.cov_scale);
    end
end