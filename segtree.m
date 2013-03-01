function [tree tree2 tree3 tree4 ucm gb_CSMG orC] = segtree(I, u, v, hu, hv, params)

M = motionMap(u, v);
rows = size(I, 1);
cols = size(I, 2);
%smooth the maps
h = fspecial('gaussian', [3 3], 0.5);
for i = 1:size(I, 3)
    I(:,:,i) = imfilter(I(:,:,i), h, 'replicate');
end

for i = 1:size(M, 3)
    M(:,:,i) = imfilter(M(:,:,i), h, 'replicate');
end

[gb_CSMG, orC] = Gb_CSMG(I, M);
ucm = contours2ucm_shugao(gb_CSMG, orC, 'imageSize'); 
params2 = params;
params2.DISPLAY = 0;
tree = segtree_from_ucm(ucm, params2); 
tree2 = cleanTree(tree, params2);
m = imfilter(sqrt(hu.*hu + hv.*hv), h, 'replicate');
tree3 = rmStatic(tree2, m, params2);
tree4 = rmRect(tree3, rows,cols, params2);

if params.DISPLAY > 0
    subplot(2, 3, 1); imagesc(gb_CSMG); colorbar;
    subplot(2, 3, 2); imagesc(ucm); colorbar;
    subplot(2, 3, 3); displayTree(I, tree);
    subplot(2, 3, 4); displayTree(I, tree2);
    subplot(2, 3, 5); displayTree(I, tree3);
    subplot(2, 3, 6); displayTree(I, tree4);
    %keyboard;
end

tree = treeDepth(tree);
tree2 = treeDepth(tree2);
tree3 = treeDepth(tree3);
tree4 = treeDepth(tree4);