function C = getRootCandidate(CC, fmap, imgs, params)
% GETROOTCANDIDATE  Generate candiate root regions from initial root
% candidates by iterative merging and non-maximum suppression

fid = [CC(:).fid];
[rows cols nfms] = size(fmap);
N = rows * cols;
C = [];


for i = 1:nfms
    C1 = CC(fid == i);
    if length(C1) < 2 
        continue;
    end
    nc = length(C1);
    T = zeros(nc);
    % init the overlap table
    for j = 1:nc-1
        for k = j+1:nc
            overlap = length(intersect(C1(j).PixelIdxList, C1(k).PixelIdxList));
            T(j, k) = overlap / min(length(C1(j).PixelIdxList), length(C1(k).PixelIdxList)); 
        end
    end
    % consequtively merge candidates to produce new candidates
    [maxovlp id] = max(T(:));
    while maxovlp >= params.er_ovlp_thre
        [y x] = ind2sub(size(T), id);
        newC.sid = [C1(y).sid C1(x).sid];
        newC.fid = i;
        newC.PixelIdxList = union(C1(y).PixelIdxList, C1(x).PixelIdxList);
        newC.Selected = 1;
        nc = length(C1);
        T(y, :) = 0;
        T(:, y) = 0;
        T(x, :) = 0;
        T(:, x) = 0;
        T1 = [T zeros(nc, 1); zeros(1, nc+1)];
        len = length(newC.PixelIdxList);
        for j = 1:nc
            if ~any(T(j,:))
                continue;
            end
            overlap = length(intersect(C1(j).PixelIdxList, newC.PixelIdxList));
            T1(j, nc+1) = overlap / min(len, length(C1(j).PixelIdxList)); 
        end
        C1 = [C1 newC];
        T = T1;
        [maxovlp id] = max(T(:));
    end
    
    % create candidate for each connected component
    bw = zeros(rows, cols);
    for j = 1:length(C1)
        bw(C1(j).PixelIdxList) = 1;
    end
    CC2 = bwconncomp(bw, 8);
    for j = 1:length(CC2)
        newC.sid = -1;
        newC.fid = i;
        newC.PixelIdxList = CC2(j).PixelIdxList;
        newC.Selected = 1;
        C1 = [C1 newC];
    end
    
    C = [C C1];
end


% perform non-maximum suppression to reduce redundant candidates
for i = 1:length(C)
    C(i).Score = sum(fmap(N * (C(i).fid - 1) + C(i).PixelIdxList)) / (length(C(i).PixelIdxList) + params.er_szprior);
end
Selected = ones(size(C));
for i = 1:nfms
    cidx = find(fid == i);
    C1 = C(cidx);
    if length(C1) < 2 
        continue;
    end
    nc = length(C1);
    T = inf(nc);
    % init the difference table
    for j = 1:nc-1
        len = length(C1(j).PixelIdxList);
        for k = j+1:nc
            len2 =  length(C1(k).PixelIdxList);
            intLen = length(intersect(C1(j).PixelIdxList, C1(k).PixelIdxList)); 
            T(j, k) = max((len - intLen)/len, (len2 - intLen)/len2);
        end
    end
    [mindiff id] = min(T(:));
    while mindiff <= params.er_diff_thre
        [y x] = ind2sub(size(T), id);          
        if C1(y).Score < C1(x).Score
            Selected(cidx(y)) = 0;
            T(y,:) = inf;
            T(:,y) = inf;
        else
            Selected(cidx(x)) = 0;
            T(x,:) = inf;
            T(:,x) = inf;
        end
        [mindiff id] = min(T(:));
    end
end

C = C(Selected > 0);

% extract feature of each candidate
for i = 1:length(C)
	im = imgs(:,:,:,C(i).fid);
    if size(C(i).PixelIdxList, 1) > 1
        C(i).PixelIdxList = C(i).PixelIdxList';
    end
	C(i).Feat = double([im(C(i).PixelIdxList); im(C(i).PixelIdxList + N); im(C(i).PixelIdxList + 2 * N)]);
end




