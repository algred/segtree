function [org_tree tree] = segtree_with_prior(I, u, v, hu, hv, fmap, params)

M = motionMap(u, v);

%smooth the maps
h = fspecial('gaussian', [3 3], 0.5);
for i = 1:size(I, 3)
    I(:,:,i) = imfilter(I(:,:,i), h, 'replicate');
end

for i = 1:size(M, 3)
    M(:,:,i) = imfilter(M(:,:,i), h, 'replicate');
end

se = strel('disk', 10, 8);
fmap = imclose(fmap, se);

[gb_CSMPG, orC] = Gb_CSMPG(I, M, fmap);
ucm = contours2ucm_shugao(gb_CSMPG, orC, 'imageSize'); 
params2 = params;
params2.DISPLAY = 0;
tree = segtree_from_ucm(ucm, params2); 
tree = saliencyTree(tree, I);
tree2 = cleanTree(tree, params2);
m = sqrt(hu.*hu + hv.*hv);
tree3 = rmStatic(tree2, m, params2);

if params.DISPLAY > 0
    h = figure('Name', 'Plots');
    set(0,'CurrentFigure',h);
    subplot(2, 3, 1); imshow(I);
    subplot(2, 3, 2); imagesc(gb_CSMPG); colorbar;
    subplot(2, 3, 3); imagesc(ucm); colorbar;
    subplot(2, 3, 4); displayTree(I, tree);
    subplot(2, 3, 5); displayTree(I, tree2);
    subplot(2, 3, 6); displayTree(I, tree3);
    pause;
end

tree = treeDepth(tree3);
org_tree = treeDepth(tree2);

tree(1).Feat = [];
tree(2:end) = treeFeat(tree(2:end), I, params);

org_tree(1).Feat = [];
org_tree(2:end) = treeFeat(org_tree(2:end), I, params);