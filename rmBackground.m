function new_trees = rmBackground(trees, fmap, params)
[rows cols nfms] = size(fmap);
N = rows * cols;
new_trees = cell(size(trees));

for ii = 1:nfms
    thisTree = trees{ii};
    PP = [thisTree(:).Parent];
    idxx = find(PP(2:end) == 1)+1; % the root is the whole frame
    dtt = treeDecendt(thisTree);
    ntt = length(thisTree);
    S = zeros(ntt, 1);
    thisOffset = (ii - 1) * N;
    for jj = 1:ntt
        S(jj) = median(fmap(thisOffset + thisTree(jj).PixelIdxList));
    end
    keep = ones(size(thisTree));
    for jj = 1:length(idxx)
        idd = idxx(jj);
        idxx2 = [idd dtt{idd}];
        if ~any(S(idxx2) > params.er_score_thre)
            keep(idxx2) = 0;
        end
    end
    % restore parent / child pointers
    idxx = find(keep);
    for jj = 2:length(idxx)
        PP2 = find(PP == idxx(jj));
        for kk = 1:length(PP2)
            thisTree(PP2(kk)).Parent = jj;
        end
    end
    new_trees{ii} = thisTree(keep > 0);
end

end