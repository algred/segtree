function vtree = treeMerge(vtree, tree, t, MC)
% TREEMERGE merges a tree at index t to a vtree. MC is the node matching 
%
% Policy for merging: add the image patch of matched nodes of tree to
% corresponding node's image patch list at index t in vtree; create new
% node in vtree for the unmatched decendents of matched nodes of tree

for i = 1:size(MC, 1)
    vtree(MC(i, 1)).list{t} = tree(MC(i, 2));
end


end