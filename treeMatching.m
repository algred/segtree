function MC = treeMatching(tree1, tree2, S)
% TREEMATCHING  This function matches two trees given the similarity matrix
% among two trees' nodes. Note tree1 and tree2 can be tree or vtree
%
% Input: tree1 and tree2 are two trees with A and B nodes respectively; S
% is a boolean matrix such that S(i, j) == 1 if node i in tree1 can be matched 
% to node j in tree2, otherwise S(i, j) == 0. S(1, :) = 0 for the special
% node of root.
%
% Output: Each row of MC records a pair of matching nodes from the two
% trees
%
% Author: Shugao Ma, 1/29/2013

S(1, :) = 0;
S(:, 1) = 0;

V = [];
for i = 2:size(S, 1)
    v = find(S(i,:));
    V = [V; [ones(length(v), 1) * i v']]; 
end

d1 = treeDecendt(tree1);
d2 = treeDecendt(tree2);

D1 = zeros(length(tree1));
D2 = zeros(length(tree2));
for i = 1:length(tree1)
    if isempty(d1{i})
        continue;
    end
    D1(i, d1{i}) = 1;
end
for i = 1:length(tree2)
    if isempty(d2{i})
        continue;
    end
    D2(i, d2{i}) = 1;
end
D1 = D1 > 0;
D2 = D2 > 0;
    

E = zeros(size(V, 1));

try
for i = 1:size(V, 1)-1
    for j = i+1:size(V, 1)
        if V(i, 1) == V(j, 1) || V(i, 2) == V(j, 2) 
            continue;
        end
        
        t1 = D1(V(j, 1), V(i, 1));
        t2 = D1(V(i, 1), V(j, 1));
        t3 = D2(V(j, 2), V(i, 2));
        t4 = D2(V(i, 2), V(j, 2));
        
        if (~t1 && t3) || (~t3 && t1) || (~t2 && t4) || (~t4 && t2)
            continue;
        end
        
        E(i, j) = 1;
        E(j, i) = 1;
    end
end
catch exception
    getReport(exception)
    keyboard;
end

Q = maximalCliques(E);
[maxv maxid] = max(sum(Q));

% Just pick the largest one. More elaborate methods can be considered in next
% version
MC = V(Q(:, maxid) > 0, :);

end