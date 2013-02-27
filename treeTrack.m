function [C tmap] = treeTrack(fmap, imgs, rootSegs, params)
% ESTIMATEROOT    Estimate the human area given the foreground estimation and candidate positions
[rows cols nfms] = size(fmap);
N = rows * cols;
C = struct;
ind = 1;
for ii = 1:nfms
    thisSegs = rootSegs{ii};
    im = imgs(:,:,:,ii);
    for jj = 1:length(thisSegs)
        C(ind).fid = ii;
        C(ind).PixelIdxList = thisSegs(jj).PixelIdxList;
        C(ind).tid = thisSegs(jj).TemplateId;
        C(ind).score = sum(fmap(N * (ii - 1) + C(ind).PixelIdxList));
        C(ind).Feat = double([im(C(ind).PixelIdxList); im(C(ind).PixelIdxList + N); im(C(ind).PixelIdxList + 2 * N)]);
        ind = ind + 1;
    end
end

T = [C(:).tid];
F = [C(:).fid];
tids = unique(T);
tmap = zeros(length(tids), 2);
tmap(:,1) = tids';
tmap(:,2) = tids';
if length(tids) <= 1
    return;
end

nt = length(tids);
L = zeros(1, nt);
for i = 1:nt
    idx = find(T == tids(i));
    L(i) = length(idx);
    F1 = F(idx);
    [mint minid] = min(F1);
    [maxt maxid] = max(F1);
    D{i} = [mint:maxt];
    I(i,:) = [idx(minid) idx(maxid)];
end

% initialize the matching table
M = zeros(nt);
for i = 1:nt-1
    mint = min(D{i});
    maxt = max(D{i});
    for j = i+1:nt
        if mint > max(D{j})
            M(i, j) = getSimilarity(C(I(i,1)).Feat, C(I(j,2)).Feat, 1:3, params.er_cov_disturb, params.er_cov_scale);
        elseif maxt < min(D{j})
            M(i, j) = getSimilarity(C(I(i,2)).Feat, C(I(j,1)).Feat, 1:3, params.er_cov_disturb, params.er_cov_scale);
        end
        M(j, i) = M(i, j);
    end
end

% merging
[L idx] = sort(L, 'descend');
for i = 1:length(L)
    tid = idx(i);
    [maxm tid2] = max(M(tid, :));
    while maxm > params.er_match_thre
        M(tid2, :) = -1;
        M(:, tid2) = -1;
        
        D{tid} = [D{tid} D{tid2}];
        T( T == tids(tid2) ) = tids(tid);
        tmap(tmap(:,1) == tids(tid2),2) = tids(tid);
        for j = 1:nt
            try
                if M(j, 1) < 0
                    continue;
                elseif j == tid || ~isempty(intersect(D{tid}, D{j}))
                    M(tid, j) = 0;
                    M(j, tid) = 0;
                    continue;
                end
                d = pdist2(D{tid}', D{j}');
                [~,ind] = min(d(:));
                [yid xid] = ind2sub(size(d), ind);
                C1 = C(T == tids(j) & F == D{j}(xid));
                C2 = C(T == tids(tid) & F == D{tid}(yid));
%                sz_ratio = length(C1.PixelIdxList) / length(C2.PixelIdxList);
%                 if sz_ratio < 1/params.er_sz_chg_thre || sz_ratio > params.er_sz_chg_thre
%                     M(tid, j) = 0;
%                 else
%                     M(tid, j) = getSimilarity(C1.Feat, C2.Feat, 1:3, params.er_cov_disturb, params.er_cov_scale);
%                 end
                M(tid, j) = getSimilarity(C1.Feat, C2.Feat, 1:3, params.er_cov_disturb, params.er_cov_scale);
                M(j, tid) = M(tid, j);
            catch exception
                getReport(exception)
                keyboard;
            end
        end
        [maxm tid2] = max(M(tid, :));
    end
end

% remove short tracks
tids = unique(T);
select = ones(size(T));
nthre = round(nfms * params.er_min_len_ratio);
L = zeros(size(tids));
for i = 1:length(tids)
    flg = (T == tids(i));
    L(i) = sum(flg);
    if L(i) < nthre
        select(flg) = 0;
    end
end
if ~any(flg)
    [~, id] = max(L);
    flg = (T == tids(id));
    select(flg) = 1;
end
select = select > 0;
C = C(select);
T = T(select);

% return result
for i = 1:length(T)
    C(i).tid = T(i);
end







