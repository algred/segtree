function decendt = treeDecendt(tree)
% TREECHILDS   Finds the indices of decendents for each tree node
%
% This assumes the childern nodes always have larger indices than parent
% nodes

P = [tree(:).Parent];

for i = length(P):-1:2
    decendt{i} = find(P == i);
    P( P == i ) = P(i);
end

decendt{1} = find(P == 1);

end