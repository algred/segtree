function displayVTree(imgs, vtree, DISPLAY_BB)

nullNode = vtree(1).list{1};
nullNode.BoundingBox = [];
nullNode.PixelIdxList = [];
nullNode.PixelIdxListBB = [];

for i = 1:size(imgs, 4)
    clear tree;
    for j = 1:length(vtree)
        if length(vtree(j).list) < i
            nullNode.Parent = vtree(j).Parent; 
            tree(j) = nullNode;
            continue;
        end
        
        if isempty(vtree(j).list{i})
            nullNode.Parent = vtree(j).Parent; 
            tree(j) = nullNode;
            continue;
        end
        
        tree(j) = vtree(j).list{i};
        tree(j).Parent = vtree(j).Parent; 
    end
    figure; displayTree2(imgs(:,:,:,i), tree, DISPLAY_BB);
end

end