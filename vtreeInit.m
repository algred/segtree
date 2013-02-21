function vtree = vtreeInit(tree, t)
% VTREEINIT  initialize a vtree from a tree

for i = 1:length(tree)
    vtree(i).list{t} = tree(i);
    vtree(i).Parent = tree(i).Parent;
end

end

