function vtree = vtreeGrow(vtree, tree, t, params)
% VTREEGROW  This function grows a vtree by merging a tree at time t to vtree
%
% Author: Shugao Ma, 1/29/2013

S = zeros(length(vtree), length(tree));

for i = 2:length(vtree)    
    for j = 2:length(tree)
    	s = zeros(size(vtree(i).list));
    	for k = 1:length(vtree(i).list)
    		if isempty(vtree(i).list{k})
    			continue;
    		end
    		s(k) = getSimilarity(tree(j), vtree(i).list{k});
    	end
        S(i, j) = max(s(k));
    end
end

S = S > params.sim_thre;

MC = treeMatching(vtree, tree, S);
vtree = treeMerge(vtree, tree, t, MC);

end