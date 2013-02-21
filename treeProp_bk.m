function [C S] = treeProp(trees, fmap, imgs, flow, motionMag, params)

%% Propagate a tree to its neighboring frames
try
    [rows cols nfms] = size(fmap);
    N = rows * cols;
    
    % remove subtrees that are on background
    trees = rmBackground(trees, fmap, params);
    
    % pre-compute motion and color histogram index maps
    [~, cmap] = rgb2ind(lab2uint8(rgb2lab(imgs(:,:,:,1))), params.er_cnbins);
    colorMaps = zeros(rows, cols, nfms);
    motionMaps = zeros(rows, cols, nfms);
    u = zeros(rows, cols, nfms);
    v = zeros(rows, cols, nfms);
    bu = zeros(rows, cols, nfms);
    bv = zeros(rows, cols, nfms);
    er_mnbins = round(params.tr_maxm/2);
    for i = 1:nfms
        colorMaps(:,:,i) = rgb2ind(lab2uint8(rgb2lab(imgs(:,:,:,i))), cmap) + 1;
        motionMaps(:,:,i) = min(er_mnbins, round(motionMag(:,:,i)/2)+1);
        u(:,:,i) = medfilt2(flow.u(:,:,i), [3 3], 'symmetric');
        v(:,:,i) = medfilt2(flow.v(:,:,i), [3 3], 'symmetric');
        bu(:,:,i) = medfilt2(flow.bu(:,:,i), [3 3], 'symmetric');
        bv(:,:,i) = medfilt2(flow.bv(:,:,i), [3 3], 'symmetric');
    end
    
    % propagate the root of each first level subtree of every frame to its
    % neighboring frames
    prm.se = strel('disk', 20, 8);
    subtrees = [];
    rootSegs = cell(1, nfms);
    for i = 1:nfms
        fmin = max(1, i - params.er_winsz);
        fmax = min(nfms, i + params.er_winsz);
        tree = trees{i};
        color_map = colorMaps(:,:,i);
        P = [tree(:).Parent];
        dt = treeDecendt(tree);
        idx = find(P(2:end) == 1);
        if isempty(idx)
            continue;
        end
        idx = idx + 1;
        nt = length(idx);
        fidx = [];
        tid = zeros(1, nt);
        for j = 1:nt
            id = idx(j);
            subtree.tree_id = i;
            subtree.node_id = id;
            subtrees = [subtrees subtree];
            tid(j) = length(subtrees);
            ndt = length(dt{id});
            rootSeg.PixelIdxList = tree(id).PixelIdxList;
            rootSeg.TemplateId = tid(j);
            rootSegs{i} = [rootSegs{i} rootSeg];
            fidx = union(fidx, tree(id).PixelIdxList);
            % color prior (keep the same during tracking)
            prm.cp = zeros(1 + ndt, params.er_cnbins);
            color_counts = histc(color_map(tree(id).PixelIdxList), 1:params.er_cnbins);
            prm.cp(1,:) = color_counts / max(color_counts);
            for k = 1:ndt
                id1 = dt{id}(k);
                color_counts = histc(color_map(tree(id1).PixelIdxList), 1:params.er_cnbins);
                prm.cp(k+1,:) = color_counts / max(color_counts);
            end
            prm.PixelIdxList = tree(id).PixelIdxList;
            % propagate forward
            for fid = i+1:fmax
                % spatial prior
                [ys xs] = ind2sub([rows, cols], prm.PixelIdxList);
                xs = round(xs + u(prm.PixelIdxList + (fid-2)*N));
                ys = round(ys + v(prm.PixelIdxList + (fid-2)*N));
                xyflg = (xs >= 1 & xs <= cols & ys >= 1 & ys <= rows);
                prm.PredPixelIdxList = sub2ind([rows, cols], ys(xyflg), xs(xyflg));
                % motion prior
                motion_counts = histc(motionMaps((fid-2) * N + prm.PixelIdxList), 1:er_mnbins);
                prm.mp = motion_counts / max(motion_counts);
                % track the segment
                r = trackSegNested;
                if isempty(r)
                    break;
                end
                rootSeg.PixelIdxList = r;
                rootSeg.TemplateId = tid(j);
                rootSegs{fid} = [rootSegs{fid} rootSeg];
            end
            % propagate backward
            for fid = i-1:-1:fmin
                % spatial prior
                [ys xs] = ind2sub([rows, cols], prm.PixelIdxList);
                xs = round(xs + bu(prm.PixelIdxList + fid*N));
                ys = round(ys + bv(prm.PixelIdxList + fid*N));
                xyflg = (xs >= 1 & xs <= cols & ys >= 1 & ys <= rows);
                prm.PredPixelIdxList = sub2ind([rows, cols], ys(xyflg), xs(xyflg));
                % motion prior
                motion_counts = histc(motionMaps((fid-2) * N + prm.PixelIdxList), 1:er_mnbins);
                prm.mp = motion_counts / max(motion_counts);
                % track the segment
                r = trackSegNested;
                if isempty(r)1:length(subtrees);
                    break;
                end
                rootSeg.PixelIdxList = r;
                rootSeg.TemplateId = tid(j);
                rootSegs{fid} = [rootSegs{fid} rootSeg];
            end
        end
        % produce root candidates from each connected component of the
        % current roots
        bw = zeros(rows, cols);
        bw(fidx) = 1;
        cc = bwconncomp(bw > 0);
        for j = 1:length(cc)
            rootSeg.PixelIdxList = cc(j).PixelIdxList;
            rootSeg.TemplateId = [];1:length(subtrees);
            for k = 1:nt1:length(subtrees);
                tmp = intersect(cc(j).PixelIdxList, tree(idx(k)).PixelIdxList);
                if isempty(tmp)
                    continue;
                end
                if length(cc(j).PixelIdxList) == length(tmp)
                    rootSeg.TemplateId = [];
                    break;
                end
                rootSeg.TemplateId = [rootSeg.TemplateId tid(j)];
            end
            if ~isempty(rootSeg.TemplateId)
                rootSegs{i} = [rootSegs{i} rootSeg];
            end
        end
    end
    
    % refine the root candidate segments
    nst = length(subtrees);
    ST = zeros(nst);
    for i = 1:nfms
        thisRootSegs = rootSegs{i};
        offset = (i-1) * N;
        for j = 1:length(thisRootSegs)
            thisRootSegs(j).Score = sum(fmap(offset + thisRootSegs(j).PixelIdxList));
        end
        rootSegs{i} = rmRedundant;
        rootSegs{i} = merge;
    end
    
    
    
catch exception
    getReport(exception)
    keyboard;
end

    function seg = trackSegNested
        % TRACKSEG   Using color, motion and spatial priors to track a segment
        %
        % fields of prm:
        %   cp, mp, color_map, motion_map, x0, y0, sp_sz, PixelIdxList
        %
        try
            % compute spatial prior
            spatial_pm = zeros(rows, cols);
            spatial_pm(prm.PredPixelIdxList) = 1;
            spatial_pm = ni(imdilate(imclose(spatial_pm, prm.se), prm.se2));
            
            % motion prior
            motion_pm = ni(medfilt2(reshape(prm.mp(motionMaps((fid-1)*N + 1: fid*N)), rows, cols), [5 5], 'symmetric'));
            % color prior
            color_pm = zeros(rows, cols);
            for ii = 1:size(prm.cp, 1)
                cp = prm.cp(ii,:);
                color_pm = color_pm + medfilt2(reshape(cp(colorMaps((fid-1)*N + 1: fid*N)), rows, cols), [5 5], 'symmetric');
            end
            color_pm = ni(color_pm);
            
            % find the match
            pm = ni(color_pm.*motion_pm.*spatial_pm);
            pm2 = round(pm * 20);cccSel
            tsz = length(prm.PixelIdxList);
            minszdiff = inf;
            for ii = 0:20
                bw = imclose((pm2 >= ii), prm.se);
                CC = bwconncomp(bw);
                for jj = 1:CC.NumObjects
                    szdiff = abs(length(CC.PixelIdxList{jj}) - tsz);
                    if szdiff < minszdiff
                        minszdiff = szdiff;
                        seg = CC.PixelIdxList{jj};
                    end
                end
            end
            
            tsz1 = length(seg);
            if tsz1 <= tsz * 0.7 || tsz1 >= tsz * 1.5
                seg = [];
                return;
            end
        catch exception1
            getReport(exception1)
            keyboard;
        end
    end


    function r = merge
        % MERGE  iterative merging root candidates
        r = thisRootSegs;
        if length(thisRootSegs) < 2
            return;
        end
        nc = length(thisRootSegs);
        T = zeros(nc);
        % init the overlap table
        for jj = 1:nc-1
            for kk = jj+1:nc
                overlap = length(intersect(thisRootSegs(jj).PixelIdxList, thisRootSegs(kk).PixelIdxList));
                T(jj, kk) = overlap / min(length(thisRootSegs(jj).PixelIdxList), length(thisRootSegs(kk).PixelIdxList));
            end
        end
        % consequtively merge candidates to produce new candidates
        [maxovlp idd] = max(T(:));
        while maxovlp >= params.er_ovlp_thre
            [y x] = ind2sub(size(T), idd);
            newC.TemplateId = [r(y).TemplateId r(x).TemplateId];
            newC.PixelIdxList = union(r(y).PixelIdxList, r(x).PixelIdxList);
            nc = length(r);
            T(y, :) = 0;
            T(:, y) = 0;rootSegs{i} = rmRedundant;
            T(x, :) = 0;
            T(:, x) = 0;
            T1 = [T zeros(nc, 1); zeros(1, nc+1)];
            len = length(newC.PixelIdxList);
            for jj = 1:nc
                if ~any(T(jj,:))
                    continue;
                end
                overlap = length(intersect(r(jj).PixelIdxList, newC.PixelIdxList));
                T1(jj, nc+1) = overlap / min(len, length(r(j1:length(subtrees);j).PixelIdxList));
            end
            r = [r newC];
            T = T1;
            [maxovlp idd] = max(T(:));
        end
    end

    function r = rmRedundant
        % MERGE  remove redundant root candidates
        r = thisRootSegs;
        if length(thisRootSegs) < 2
            return;
        end
        nc = length(thisRootSegs);
        T = inf(nc);
        % compute the difference table
        for jj = 1:nc-1
            len = length(r(jj).PixelIdxList);
            for kk = jj+1:nc
                len2 = length(r(kk).PixelIdxList);
                intLen = length(intersect(r(jj).PixelIdxList, r(kk).PixelIdxList)); 
                T(jj, kk) = max((len - intLen)/len, (len2 - intLen)/len2);
            end
        end
        T = (T <= params.er_diff_thre);
        if ~any(T(:))
            return;
        end
        ccc = bwconncomp(T);
        rootSel = ones(size(thisRootSegs));
        thisScore = [thisRootSegs(:).Score];
        for jj = 1:length(ccc)  
            [yys xxs] = ind2sub([nc nc], ccc(jj).PixelIdxList);
            cccIdx = unique([yys xxs]);
            tpid = [thisRootSegs(cccIdx).TemplateId];
            tpid = sort(tpid);
            for kk = 1:length(tpid)-1
                for hh = 2:length(tpid)
                    ST(kk, hh) = 1;
                end
            end
            cccScore = thisScore(cccIdx);
            [~, bestccc] = max(cccScore);
            rootSel(cccIdx([1:bestccc-1 bestccc+1:end])) = 0;
        end
    end
end