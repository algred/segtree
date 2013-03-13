function [rootSegs subtrees]  = treeProp2(trees, fmap, imgs, colorMaps, colorPriors, u, v, bu, bv,  params)

%% Propagate a tree to its neighboring frames
try
    [rows cols nfms] = size(fmap);
    N = rows * cols;
    min_root_sz = round(params.min_root_sz_ratio * N);
    
    % propagate the root of each first level subtree of every frame to its
    % neighboring frames
    subtrees = [];
    rootSegs = cell(1, nfms);
    er_grid_r = min(cols / 50, rows / 50);
    for i = 1:nfms
        if isempty(trees{i})
            continue;
        end
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
        tid = zeros(1, nt);
        if params.verbose > 0
            fprintf('Processing frame %d, subtree number = %d\n', i, nt);
        end
        for j = 1:nt
            id = idx(j);
            if length(tree(id).PixelIdxList) < min_root_sz
                continue;
            end
            subtree.tree_id = i;
            subtree.node_id = id;
            tid(j) = length(subtrees)+1;
            ndt = length(dt{id});
            if size(tree(id).PixelIdxList, 1) > 1
                tree(id).PixelIdxList = tree(id).PixelIdxList';
            end
            rootSeg.FromProp = 0;
            rootSeg.PixelIdxList = tree(id).PixelIdxList;
            rootSeg.TemplateId = tid(j);
            rootSegs{i} = [rootSegs{i} rootSeg];
            
            % color prior (keep the same during tracking)
            prm.cp = zeros(1 + ndt, params.er_cnbins);
            color_counts = histc(color_map(tree(id).PixelIdxList), 1:params.er_cnbins);
            prm.cp(1,:) = color_counts / sum(color_counts);
            for k = 1:ndt
                id1 = dt{id}(k);
                color_counts = histc(color_map(tree(id1).PixelIdxList), 1:params.er_cnbins);
                prm.cp(k+1,:) = color_counts / sum(color_counts);
            end
            [ys xs] = ind2sub([rows, cols], tree(id).PixelIdxList);
            bbox = [i min(xs) min(ys) max(xs)-min(xs)+1 max(ys)-min(ys)+1];
            if params.verbose > 0
                fprintf('propagate subtree %d to frame: ', j);
            end
            % propagate forward
            prm.PixelIdxList = tree(id).PixelIdxList;
            for fid = i+1:fmax
                tsz = length(prm.PixelIdxList);
                % spatial prior
                [ys xs] = ind2sub([rows, cols], prm.PixelIdxList); 
                xs = round(xs + u(prm.PixelIdxList + (fid-2)*N));
                ys = round(ys + v(prm.PixelIdxList + (fid-2)*N));
                [bb prm.SpatialPM] = idx2predMask;
                if isempty(bb)
                    break;
                end                
                % track the segment
                r = trackSegNested;
                if isempty(r)
                    break;
                end
                if params.DISPLAY > 0
                    im = imgs(:,:,:,fid);
                    im2 = imgs(:,:,:,fid-1);
                    mask1 = uint8(zeros(rows, cols,3));
                    mask2 = uint8(zeros(rows, cols,3));
                    mask1(r) = im(r); mask1(r+N) = im(r+N); mask1(r+2*N) = im(r+2*N);
                    mask2(r) = im2(r); mask2(r+N) = im2(r+N); mask2(r+2*N) = im2(r+2*N);
                    subplot(2, 1, 1); imshow(mask2); subplot(2, 1, 2); imshow(mask1);
                    keyboard;
                end
                rootSeg.FromProp = 1;
                rootSeg.PixelIdxList = r;
                rootSeg.TemplateId = tid(j);
                rootSegs{fid} = [rootSegs{fid} rootSeg];
                prm.PixelIdxList = r;
                [ys xs] = ind2sub([rows, cols], prm.PixelIdxList);
                bbox = [fid min(xs) min(ys) max(xs) - min(xs) + 1  max(ys) - min(ys) + 1; bbox];
                if params.verbose > 0
                    fprintf('%d ', fid);
                end
            end
            % propagate backward
            prm.PixelIdxList = tree(id).PixelIdxList;
            for fid = i-1:-1:fmin
                tsz = length(prm.PixelIdxList);
                % spatial prior
                [ys xs] = ind2sub([rows, cols], prm.PixelIdxList);
                xs = round(xs + bu(prm.PixelIdxList + fid*N));
                ys = round(ys + bv(prm.PixelIdxList + fid*N));
                [bb prm.SpatialPM] = idx2predMask;
                if isempty(bb)
                    break;
                end
                % track the segment
                r = trackSegNested;
                if isempty(r)
                    break;
                end
                if params.DISPLAY > 0
                    im = imgs(:,:,:,fid);
                    im2 = imgs(:,:,:,fid+1);
                    mask1 = uint8(zeros(rows, cols,3));
                    mask2 = uint8(zeros(rows, cols,3));
                    mask1(r) = im(r); mask1(r+N) = im(r+N); mask1(r+2*N) = im(r+2*N);
                    mask2(r) = im2(r); mask2(r+N) = im2(r+N); mask2(r+2*N) = im2(r+2*N);
                    subplot(2, 1, 1); imshow(mask2); subplot(2, 1, 2); imshow(mask1);
                    keyboard;
                end
                rootSeg.FromProp = 1;
                rootSeg.PixelIdxList = r;
                rootSeg.TemplateId = tid(j);
                rootSegs{fid} = [rootSegs{fid} rootSeg];
                prm.PixelIdxList = r;
                [ys xs] = ind2sub([rows, cols], prm.PixelIdxList);
                bbox = [fid min(xs) min(ys) max(xs) - min(xs) + 1  max(ys) - min(ys) + 1; bbox];
                if params.verbose > 0
                    fprintf('%d ', fid);
                end
            end
            subtree.bbox = bbox;
            subtrees = [subtrees subtree];
            if params.verbose > 0
               fprintf('\n');
            end
        end
    end
    
    % refine the root candidate segments
    nst = length(subtrees);
    sfid = [subtrees(:).tree_id];
    CF = zeros(nst); % conflict matrix
    for i = 1:nfms
        flg = (sfid == i);
        if ~any(flg)
            continue;
        end
        CF(flg, flg) = 1;
    end
    CF(sub2ind(size(CF), [1:nst], [1:nst])) = 0; 
    
    ST = eye(nst);
    for i = 1:nfms
        thisRootSegs = rootSegs{i};
        offset = (i-1) * N;
        for j = 1:length(thisRootSegs)
            thisRootSegs(j).Score = sum(fmap(offset + thisRootSegs(j).PixelIdxList)) / (length(thisRootSegs(j).PixelIdxList)+params.er_szbias);
        end
        mergeTracks;
        rootSegs{i} = thisRootSegs;
    end
    
    % identify tracks
    ST = double(ST > 0);
    trackID = connComp(ST);
    for i = 1:length(subtrees)
        subtrees(i).track_id = trackID(i);
    end

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
        fromPropFlg = [rootSegs{i}(:).FromProp];
        for j = 1:length(rtglabel)
            rtidx = find(rtgidx == rtglabel(j));
            fpflg = fromPropFlg(rtidx) < 1;
            if ~any(fpflg)
                [~, bestrt] = max(rtscore(rtidx));
                rtsel(rtidx([1:bestrt-1 bestrt+1:end])) = 0;
            else
                rtsel(rtidx(~fpflg)) = 0;
                if(sum(fpflg) > 1)
                    keyboard;
                end
            end
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

    function seg = trackSegNested
        % TRACKSEG   Using color, motion and spatial priors to track a segment
        %
        % fields of prm:
        %   cp, mp, color_map, motion_map, x0, y0, sp_sz, PixelIdxList
        %
        try
            % color prior
            color_pm = zeros(bb(4), bb(3));
            for ii = 1:size(prm.cp, 1)
                cp = prm.cp(ii,:)./colorPriors(fid, :);
                thisColorMap = colorMaps(bb(2):bb(2)+bb(4)-1, bb(1):bb(1)+bb(3)-1, fid);
                color_pm = color_pm + reshape(cp(thisColorMap(:)), bb(4), bb(3));
            end
            color_pm = ni(color_pm);
            if params.DISPLAY > 0
                subplot(2, 1, 1); imagesc(color_pm);colorbar;
                subplot(2, 1, 2); imagesc(prm.SpatialPM);colorbar;
                keyboard;
            end
            % find the match
            pm = ni(color_pm.*prm.SpatialPM);
            pm2 = round(pm * 20);
            minszdiff = inf;
            seg = [];
            for ii = 1:19
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
            if tsz1 <= tsz * params.er_max_shrink || tsz1 >= tsz * params.er_max_expand
                seg = [];
                return;
            end
            [yys xxs] = ind2sub([bb(4), bb(3)], seg);
            seg = sub2ind([rows, cols], yys + bb(2) - 1, xxs + bb(1) -1);
        catch exception1
            getReport(exception1)
            keyboard;
        end
    end

    function mergeTracks
        if length(thisRootSegs) < 2
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
        gidx = connComp(T);
        glabel = unique(gidx);
        for jj = 1:length(glabel)
            tpid = [thisRootSegs(gidx == glabel(jj)).TemplateId];
            tpid = sort(tpid);
            for kk = 1:length(tpid)-1
                for hh = kk+1:length(tpid)
                    ST1 = ST;
                    ST1(tpid(kk), tpid(hh)) = 1;
                    ST1(tpid(hh), tpid(kk)) = 1;
                    if ~isConflict(ST1)
                        ST = ST1;
                    end
                end
            end
        end
    end

    function flg = isConflict(thisST)
        tracks = connComp(thisST);
        trackIDS = unique(tracks);
        flg = 0;
        for ii = 1:length(trackIDS)
            ids = find(tracks == trackIDS(ii));
            CF1 = CF(ids, ids);
            if any(CF1(:))
                flg = 1;
                break;
            end
        end
    end

    function [thisBB r] = idx2predMask
        xyflg = xs >= 1 & xs <= cols & ys >= 1 & ys <= rows;
        xs = xs(xyflg);
        ys = ys(xyflg);
        if length(xs) < 10
            r = [];
            thisBB = [];
            return;
        end
        prm.PixelIdxList = sub2ind([rows, cols], ys, xs);
        
        thisBB = [min(xs) min(ys) max(xs)-min(xs)+1 max(ys)-min(ys)+1];
        grid_width = ceil(thisBB(3) / er_grid_r);
        grid_height = ceil(thisBB(4) / er_grid_r);
        gxs = min(grid_width, floor((xs-thisBB(1)+1)/er_grid_r) + 1);
        gys = min(grid_height, floor((ys-thisBB(2)+1)/er_grid_r) + 1);
        grid = zeros(grid_height, grid_width);
        grid(sub2ind([grid_height, grid_width], gys, gxs)) = 1;
        grid = imresize(grid, [thisBB(4), thisBB(3)], 'nearest');
        [yys xxs] = find(grid);
        predBB = boundingBox(xxs', yys');
        r = double(poly2mask(predBB(1, [1:end 1]), predBB(2, [1:end 1]), thisBB(4), thisBB(3)));
        r(sub2ind([thisBB(4), thisBB(3)], yys, xxs)) = 2;
        if sum(r(:)) < (tsz*params.er_max_shrink)
            r = [];
            thisBB = [];
            return;
        end
    end

end