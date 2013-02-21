function [C S] = propCandSegs(trees, fmap, imgs, flow, motionMag, params)

%% Generate candidates in each frame
try
    [rows cols nfms] = size(fmap);
    N = rows * cols;
    
    % record all the segments in the video
    allS = [];
    for i = 1:nfms
        P = [trees{i}(:).Parent];
        nSeg = length(trees{i});
        if nSeg <= 1
            continue;
        end
        thisS = struct;
        thisScore = zeros(nSeg, 1);
        for j = nSeg:-1:2  % the root node is the whole frame
            thisS(j).fid = i;
            thisS(j).tsid = j;
            bbox = trees{i}(j).BoundingBox;
            bbw1 = sqrt((bbox(1,1) - bbox(1,2))*(bbox(1,1) - bbox(1,2)) + (bbox(2,1) - bbox(2,2))*(bbox(2,1) - bbox(2,2)));
            bbw2 = sqrt((bbox(1,2) - bbox(1,3))*(bbox(1,2) - bbox(1,3)) + (bbox(2,2) - bbox(2,3))*(bbox(2,2) - bbox(2,3)));
            if(max(bbw1, bbw2)/ min(bbw1, bbw2) > 10)
                thisScore(j) = 0;
            else
                thisScore(j) = sum(thisScore(P == j)) + sum(fmap(trees{i}(j).PixelIdxList + N*(i-1)))/length(trees{i}(j).PixelIdxList);
            end
            thisS(j).score = thisScore(j);
            thisS(j).propLeftFid = -1;
            thisS(j).canPropLeft = 1;
            thisS(j).canPropRight = 1;
            thisS(j).propRightFid = -1;
            thisS(j).propLeftCid = -1;
            thisS(j).propRightCid = -1;
            thisS(j).cp = [];
        end
        allS = [allS; thisS(2:end)'];
    end
    
    [~, idx] = sort([allS(:).score], 'descend');
    S = allS(idx);
    
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
    
    % track top segments within each sliding window
    prm.se = strel('disk', 10, 8);
    prm.se2 = strel('disk', 10, 8);
    %prm.gh = fspecial('gaussian', 10, 10);
    fidx = [S(:).fid];
    ind = 1;
    for i = 1:nfms
        fmin = max(1, i - params.er_winsz);
        fmax = min(nfms, i + params.er_winsz);
        sidx = find(fidx >= fmin & fidx <= fmax);
        if length(sidx) > params.er_k
            sidx = sidx(1:params.er_k);
        end
        % propergate every top k segment within the time window onto frame i
        for j = 1:length(sidx);
            id = sidx(j);
            seg = trees{S(id).fid}(S(id).tsid);
            decendt = treeDecendt(trees{S(id).fid});
            if S(id).propLeftFid < 0
                thisDecendt = decendt{S(id).tsid};
                C(ind).sid = id;
                C(ind).fid = S(id).fid;
                C(ind).PixelIdxList = seg.PixelIdxList;
                if size(C(ind).PixelIdxList, 1) > 1
                    C(ind).PixelIdxList = C(ind).PixelIdxList';
                end
                C(ind).Selected = 0;
                S(id).cp = zeros(1 + length(thisDecendt), params.er_cnbins);
                color_map = colorMaps(:,:,S(id).fid);
                color_counts = histc(color_map(seg.PixelIdxList), 1:params.er_cnbins);
                S(id).cp(1,:) = color_counts / max(color_counts);
                for k = 1:length(thisDecendt)
                    id1 = thisDecendt(k);
                    color_counts = histc(color_map(trees{S(id).fid}(id1).PixelIdxList), 1:params.er_cnbins);
                    S(id).cp(k+1,:) = color_counts / max(color_counts);
                end
                S(id).propLeftFid = S(id).fid;
                S(id).propRightFid = S(id).fid;
                S(id).propLeftCid = ind;
                S(id).propRightCid = ind;
                ind = ind + 1;
            end
            
            if i > S(id).propRightFid   % need to propergate forward
                if S(id).canPropRight < 1 % can't be propergated
                    continue;
                end
                prm.cp = S(id).cp;
                prm.PixelIdxList = C(S(id).propRightCid).PixelIdxList;
                for k = S(id).propRightFid + 1:i
                    % for spatial prior
                    [ys xs] = ind2sub([rows, cols], prm.PixelIdxList);
                    xs = round(xs + u(prm.PixelIdxList + (k-2)*N));
                    ys = round(ys + v(prm.PixelIdxList + (k-2)*N));
                    xyflg = (xs >= 1 & xs <= cols & ys >= 1 & ys <= rows);
                    prm.PredPixelIdxList = sub2ind([rows, cols], ys(xyflg), xs(xyflg));
                    % motion prior
                    motion_counts = histc(motionMaps((k-2) * N + prm.PixelIdxList), 1:er_mnbins);
                    prm.mp = motion_counts / max(motion_counts);
                    % track the segment
                    r = trackSegNested;
                    if isempty(r)
                        S(id).canPropRight = 0;
                        break;
                    end
                    C(ind).sid = id;
                    C(ind).fid = k;
                    C(ind).PixelIdxList = r;
                    if size(C(ind).PixelIdxList, 1) > 1
                        C(ind).PixelIdxList = C(ind).PixelIdxList';
                    end
                    C(ind).Selected = 0;
                    prm.PixelIdxList = C(ind).PixelIdxList;
                    S(id).propRightFid = k;
                    S(id).propRightCid = ind;
                    ind = ind + 1;
                end
            elseif i < S(id).propLeftFid    % need to propergate backward
                if S(id).canPropLeft < 1  % can't be propergated backward
                    continue;
                end
                prm.cp = S(id).cp;
                prm.PixelIdxList = C(S(id).propLeftCid).PixelIdxList;
                for k = S(id).propLeftFid-1:-1:i
                    % for spatial prior
                    [ys xs] = ind2sub([rows, cols], prm.PixelIdxList);
                    xs = round(xs + bu(prm.PixelIdxList + k*N));
                    ys = round(ys + bv(prm.PixelIdxList + k*N));
                    xyflg = (xs >= 1 & xs <= cols & ys >= 1 & ys <= rows);
                    prm.PredPixelIdxList = sub2ind([rows, cols], ys(xyflg), xs(xyflg));
                    % motion prior
                    motion_counts = histc(motionMaps(k*N + prm.PixelIdxList), 1:er_mnbins);
                    prm.mp = motion_counts / max(motion_counts);
                    % track the segment
                    r = trackSegNested;
                    if isempty(r)
                        S(id).canPropLeft = 0;
                        break;
                    end
                    C(ind).sid = id;
                    C(ind).fid = k;
                    C(ind).PixelIdxList = r;
                    if size(C(ind).PixelIdxList, 1) > 1
                        C(ind).PixelIdxList = C(ind).PixelIdxList';
                    end
                    C(ind).Selected = 0;
                    prm.PixelIdxList = C(ind).PixelIdxList;
                    S(id).propLeftFid = k;
                    S(id).propLeftCid = ind;
                    ind = ind + 1;
                end
            end
            
            cflg = ([C(:).sid] == id) & ([C(:).fid] == i);
            if ~any(cflg)
                continue;
            end
            C(cflg).Selected = 1;
        end
    end
    
    flg = [C(:).Selected] > 0;
    C = C(flg);
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
            motion_pm = ni(medfilt2(reshape(prm.mp(motionMaps((k-1)*N + 1: k*N)), rows, cols), [5 5], 'symmetric'));
            % color prior
            color_pm = zeros(rows, cols);
            for ii = 1:size(prm.cp, 1)
                cp = prm.cp(ii,:);
                color_pm = color_pm + medfilt2(reshape(cp(colorMaps((k-1)*N + 1: k*N)), rows, cols), [5 5], 'symmetric');
            end
            color_pm = ni(color_pm);
            
            % find the match
            pm = ni(color_pm.*motion_pm.*spatial_pm);
            pm2 = round(pm * 20);
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
end