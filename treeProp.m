function [rootSegs trackID]  = treeProp(trees, fmap, imgs, flow, motionMag, params)

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
        colorMaps(:,:,i) = medfilt2(rgb2ind(lab2uint8(rgb2lab(imgs(:,:,:,i))), cmap) + 1, [5 5], 'symmetric');
        motionMaps(:,:,i) = medfilt2(min(er_mnbins, round(motionMag(:,:,i)/2))+1, [5 5], 'symmetric');
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
    er_grid_r = min(cols / 50, rows / 50);
    grid_width = round(cols / er_grid_r);
    grid_height = round(rows / er_grid_r);
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
                prm.SpatialPM = idx2predMask;
                % motion prior
                motionX = motionMaps((fid-2) * N + prm.PixelIdxList);
                prm.MotionPM = ni(reshape(normpdf(motionMaps((fid-1) * N + 1 : fid*N), ...
                    mean(motionX), std(motionX)), rows, cols));
                % track the segment
                r = trackSegNested;
                if isempty(r)
                    break;
                end
                rootSeg.PixelIdxList = r;
                rootSeg.TemplateId = tid(j);
                rootSegs{fid} = [rootSegs{fid} rootSeg];
                prm.PixelIdxList = r;
            end
            % propagate backward
            for fid = i-1:-1:fmin
                % spatial prior
                [ys xs] = ind2sub([rows, cols], prm.PixelIdxList);
                xs = round(xs + bu(prm.PixelIdxList + fid*N));
                ys = round(ys + bv(prm.PixelIdxList + fid*N));
                prm.SpatialPM = idx2predMask;
                % motion prior
                motionX = motionMaps(fid * N + prm.PixelIdxList);
                prm.MotionPM = ni(reshape(normpdf(motionMaps((fid-1) * N + 1 : fid*N), ...
                    mean(motionX), std(motionX)), rows, cols));
                % track the segment
                r = trackSegNested;
                if isempty(r)
                    break;
                end
                rootSeg.PixelIdxList = r;
                rootSeg.TemplateId = tid(j);
                rootSegs{fid} = [rootSegs{fid} rootSeg];
                prm.PixelIdxList = r;
            end
        end
    end
    
    % refine the root candidate segments
    nst = length(subtrees);
    ST = eye(nst);
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
    ST = double(ST > 0);
    trackID = connComp(ST);

    % remove redundant again
    for i = 1:nfms
        nc = length(rootSegs{i});
        for j = 1:nc
            rootSegs{i}(j).TemplateId = trackID(rootSegs{i}(j).TemplateId);
        end
        RT = zeros(nc);
        for j = 1:nc-1
            for k = j+1:nc
                if rootSegs{i}(j).TemplateId == rootSegs{i}(k).TemplateId
                    RT(j, k) = 1;
                    RT(k, j) = 1;
                end
            end
        end
        RT = double(RT > 0);
        if ~any(RT(:))
            continue;
        end
        rtgidx = connComp(RT); 
        rtglabel = unique(rtgidx);
        rtsel = ones(1, nc);
        rtscore = [rootSegs{i}(:).Score];
        for j = 1:length(rtglabel)
            rtidx = find(rtgidx == rtglabel(j));
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
            keep = ones(size(thisTree));
            for jj = 1:length(idxx)
                idd = idxx(jj);
                idxx2 = [idd dtt{idd}];
                if ~any(S(idxx2) > params.er_score_thre)
                    keep(idxx2) = 0;
                end
            end
            % restore parent / child pointers 
            idxx = find(keep);
            for jj = 2:length(idxx)
                PP2 = find(PP == idxx(jj));
                for kk = 1:length(PP2)
                    thisTree(PP2(kk)).Parent = jj;
                end
            end
            new_trees{ii} = thisTree(keep > 0);
        end
    end

    function seg = trackSegNested
        % TRACKSEG   Using color, motion and spatial priors to track a segment
        %
        % fields of prm:
        %   cp, mp, color_map, motion_map, x0, y0, sp_sz, PixelIdxList
        %
        try
            % color prior
            color_pm = zeros(rows, cols);
            for ii = 1:size(prm.cp, 1)
                cp = prm.cp(ii,:);
                color_pm = color_pm + reshape(cp(colorMaps((fid-1)*N + 1: fid*N)), rows, cols);
            end
            color_pm = ni(color_pm);
            tsz = length(prm.PixelIdxList);
            % find the match
            pm = ni(color_pm.*prm.MotionPM.*prm.SpatialPM);
            pm2 = round(pm * 10);
            minszdiff = inf;
            seg = [];
            for ii = 2:19
                bw = imfill((pm2 >= ii), 'holes');
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
            if tsz1 <= tsz * 0.7 || tsz1 >= tsz * 1.3
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
        if length(thisRootSegs) < 2
            r = thisRootSegs;
            return;
        end
        nnc = length(thisRootSegs);
        T = inf(nnc);
        % compute the difference table
        for jj = 1:nnc-1
            len = length(thisRootSegs(jj).PixelIdxList);
            for kk = jj+1:nnc
                len2 = length(thisRootSegs(kk).PixelIdxList);
                intLen = length(intersect(thisRootSegs(jj).PixelIdxList, thisRootSegs(kk).PixelIdxList));
                T(jj, kk) = max((len - intLen)/len, (len2 - intLen)/len2);
                T(kk, jj) = T(jj, kk);
            end
        end
        T = double(T <= params.er_diff_thre);
        if ~any(T(:))
            r = thisRootSegs;
            return;
        end
        rootSel = ones(size(thisRootSegs));
        gidx = connComp(T);
        glabel = unique(gidx);
        thisScore = [thisRootSegs(:).Score];
        for jj = 1:length(glabel)
            cccIdx = find(gidx == glabel(jj));
            tpid = [thisRootSegs(cccIdx).TemplateId];
            tpid = sort(tpid);
            for kk = 1:length(tpid)-1
                for hh = kk+1:length(tpid)
                    ST(tpid(kk), tpid(hh)) = 1;
                    ST(tpid(hh), tpid(kk)) = 1;
                end
            end
            cccScore = thisScore(cccIdx);
            [~, bestccc] = max(cccScore);
            rootSel(cccIdx([1:bestccc-1 bestccc+1:end])) = 0;
        end
        r = thisRootSegs(rootSel > 0);
    end

    function r = idx2predMask
        xyflg = xs >= 1 & xs <= cols & ys >= 1 & ys <= rows;
        xs = xs(xyflg);
        ys = ys(xyflg);
        if length(xs) < 100
            r = (zeros(rows, cols) > 0);
            return;
        end
        gxs = min(grid_width, floor(xs/er_grid_r) + 1);
        gys = min(grid_height, floor(ys/er_grid_r) + 1);
        grid = zeros(grid_height, grid_width);
        grid(gys, gxs) = 1;
        grid = imresize(grid, [rows, cols], 'nearest');
        [yys xxs] = find(grid);
        predBB = boundingBox(xxs', yys');
        r = poly2mask(predBB(1, [1:end 1]), predBB(2, [1:end 1]), rows, cols);
    end
end