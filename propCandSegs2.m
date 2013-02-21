function rootSegs = propCandSegs2(trees, fmap, imgs, flow, motionMag, params)

%% Propagate a tree to its neighboring frames

try
    [rows cols nfms] = size(fmap);
    N = rows * cols;
    
    % remove subtrees that are on background
    if params.verbose > 0
        fprintf('Removing background subtrees\n');
    end
    trees = rmBackgroundNested;

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
    prm.se = strel('disk', 10, 8);
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
        if params.verbose > 0
            fprintf('Processing frame %d, subtree number = %d\n', i, nt);
        end
        for j = 1:nt
            if params.verbose > 0
                fprintf('propagate subtree %d ...\n', j);
            end
            id = idx(j);
            subtree.tree_id = i;
            subtree.node_id = id;
            subtrees = [subtrees subtree];
            tid(j) = length(subtrees);
            ndt = length(dt{id});
            if size(tree(id).PixelIdxList, 1) > 1
                tree(id).PixelIdxList = tree(id).PixelIdxList';
            end
            rootSeg.PixelIdxList = tree(id).PixelIdxList;
            rootSeg.TemplateId = tid(j);
            rootSegs{i} = [rootSegs{i} rootSeg];
            fidx = union(fidx, tree(id).PixelIdxList);
            % color prior treeProp.m(keep the same during tracking)
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
                motion_counts = histc(motionMaps(fid * N + prm.PixelIdxList), 1:er_mnbins);
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
        if params.verbose > 0
            fprintf('remove redundant on frame %d, before = %d, after = %d\n', i, length(thisRootSegs), length(rootSegs{i}));
        end
    end
    
    % identify tracks
    ST = (ST > 0);
    stcc = bwconncomp(ST);
    trackID = 1:nst;
    for i = 1:stcc.NumObjects
        [ys xs] = ind2sub([nst nst], stcc.PixelIdxList{i});
        stid = unique([ys; xs]);
        trackID(stid(2:end)) = stid(1);
    end
    % remove redundant again
    for i = 1:nfms
        nc = length(rootSegs{i});
        for j = 1:nc
            rootSegs{i}(j).TemplateId = trackID(rootSegs{i}(j).TemplateId);
        end
        RT = zeros(nc);
        for j = 1:nc-1
            for k = 2:nc
                if rootSegs{i}(j).TemplateId == rootSegs{i}(k).TemplateId
                    RT(j, k) = 1;
                end
            end
        end
        RT = RT > 0;
        if ~any(RT(:))
            continue;
        end
        rtcc = bwconncomp(RT);
        rtscore = [rootSegs{i}(:).Score];
        rtsel = ones(1, nc);
        for j = 1:rtcc.NumObjects
            [ys xs] = ind2sub([nc nc], rtcc.PixelIdxList{j});
            rtidx = unique([ys; xs]);
            [~, bestrt] = max(rtscore(rtidx));
            rtsel(rtidx([1:bestrt-1 bestrt+1:end])) = 0;
        end
        if params.verbose > 0
            fprintf('remove redundant on frame %d again, before = %d', i, length(rootSegs{i}));
        end
        rootSegs{i} = rootSegs{i}(rtsel > 0);
        if params.verbose > 0
            fprintf('after = %d \n', length(rootSegs{i}));
        end
        % produce root candidates from each connected component of the
        % current roots
        bw = zeros(rows, cols);
        bw(fidx) = 1;
        cc = bwconncomp(bw > 0);
        for j = 1:length(cc.PixelIdxList)
            rootSeg.PixelIdxList = cc.PixelIdxList{j};
            rootSeg.TemplateId = [];
            nt = length(rootSegs{i});
            for k = 1:nt;
                tmp = intersect(cc.PixelIdxList{j}, rootSeg{i}(k).PixelIdxList);
                if isempty(tmp)
                    continue;
                end
                if length(cc.PixelIdxList{j}) == length(tmp)
                    rootSeg.TemplateId = [];
                    break;
                end
                rootSeg.TemplateId = [rootSeg.TemplateId rootSeg{i}(k).TemplateId];
            end
            if ~isempty(rootSeg.TemplateId)
                rootSegs{i} = [rootSegs{i} rootSeg];
            end
        end
    end
    
catch exception
    getReport(exception)
    keyboard;
end
    function new_trees = rmBackgroundNested
        for ii = 1:nfms
            thisTree = trees{ii};
            PP = [thisTree(:).Parent];
            idxx = find(PP(2:end) == 1)+1; % the root is the whole frame
            dtt = treeDecendt(thisTree);
            ntt = length(thisTree);
            S = zeros(ntt, 1);
            thisOffset = (ii - 1) * N;
            for jj = 1:ntt
                S(jj) = sum(fmap(thisOffset + thisTree(jj).PixelIdxList))/ length(thisTree(jj).PixelIdxList);
            end
            tree1 = thisTree(1);
            for jj = 1:length(idxx)
                idd = idxx(jj);
                idxx2 = [idd dtt{idd}];
                if any(S(idxx2) > params.er_score_thre)
                    tree1 = [tree1 thisTree(idxx2)];
                end
            end
            new_trees{ii} = tree1;
        end
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
            spatial_pm = ni(imdilate(imclose(spatial_pm, prm.se), prm.se));
            
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
            pm2 = round(pm * 10);
            tsz = length(prm.PixelIdxList);
            minszdiff = inf;
            for ii = 0:10
                bw = imclose((pm2 >= ii), prm.se);
                CC = bwconncomp(bw);
                for jj = 1:CC.NumObjects
                    szdiff = abs(length(CC.PixelIdxList{jj}) - tsz);
                    if szdiff < minszdiff
                        minszdiff = szdiff;
                        seg = CC.PixelIdxList{jj}';
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



    function r = rmRedundant
        % MERGE  remove redundant root candidates
        r = thisRootSegs;
        if length(thisRootSegs) < 2
            return;
        end
        nnc = length(thisRootSegs);
        T = inf(nnc);
        % compute the difference table
        for jj = 1:nnc-1
            len = length(r(jj).PixelIdxList);
            for kk = jj+1:nnc
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
        for jj = 1:ccc.NumObjects
            [yys xxs] = ind2sub([nnc nnc], ccc.PixelIdxList{jj});
            cccIdx = unique([yys; xxs]);
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
        r = r(rootSel > 0);
    end
end