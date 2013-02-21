function tree = saliencyTree(tree, img)
% SALIENCYTREE Given a tree, compute a saliency score for each tree node

sm = saliency_map(img);

for i = 1:length(tree)
    tree(i).Saliency = sum(sm(tree(i).PixelIdxListBB))./length(tree(i).PixelIdxListBB);
end