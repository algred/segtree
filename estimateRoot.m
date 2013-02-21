function R = estimateRoot(C, params)
% ESTIMATEROOT    Estimate the human area given the foreground estimation and 
% candidate positions using dynamic programming


% Compute the dynamic programming table
T = cell(nfms, 1);
B = cell(nfms, 1);
I = cell(nfms, 1);
fid = [C(:).fid];
start = 1;
idx = find(fid == 1);
while isempty(idx)
	start = start + 1;
	idx = find(fid == start);
end
	
C1 = C(idx);
t = [C1(:).score];
%T{start} = log(t ./ sum(t));
T{start} = log(t);
B{start} = 1:length(C1);
I{start} = idx;
fidx(1) = start;
for i = start+1:nfms
	idx = find(fid == i);
	if isempty(idx)
		continue;
	end
	I{i} = idx;
	C1 = C(idx);
	C2 = C(I{fidx(end)});
	t1 = [C1(:).score];
	%t1 = t1 ./ sum(t1);
    t1 = log(t1);
	t2 = zeros(size(C1));
    b = zeros(size(C1));
	for j = 1:length(C1)
		maxscore = -inf;
		id = -1;
		for k = 1:length(C2)
			sim = log(getSimilarity(C1(j).Feat,C2(k).Feat, 1:3, params.er_cov_disturb, params.er_cov_scale));
            if sim < params.er_sim_thre
                sim = -inf;
            end
			score = T{fidx(end)}(k) + sim;
			if score > maxscore
				maxscore = score;
				id = k;
			end
        end
        if maxscore <= -inf
            t2(j) = t1(j);
            b(j) = -1;
        else
            t2(j) = t1(j) + maxscore;
            b(j) = id;
        end
	end
	T{i} = t2;
	B{i} = b;
	fidx = [fidx i];
end

[maxscore id] = max(T{fidx(end)});
nvalfms = length(fidx);
cidx = zeros(1, nvalfms);
r(fidx(end)) = id;
cidx(nvalfms) = I{fidx(end)}(id);
for i = length(fidx)-1:-1:1
	r(fidx(i)) = B{fidx(i+1)}(r(fidx(i+1)));
    if r(fidx(i)) < 0
        cidx(i) = -1;
    else
        cidx(i) = I{fidx(i)}(r(fidx(i)));
    end
end

R = C(cidx(cidx > 0));

end