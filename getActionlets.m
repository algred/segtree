function [A T] = getActionlets(trees, imgs, flow, fmap, params)
% Extract actionlets from a video using segment trees
%
% A is a struct array containing segments of actionlets with the fields:
%   tid: track id
%   pid: parent actionlet id
%   isRoot: is 1 if is from root segment
%   start: starting frame
%   end:   ending frame
%   bbox: horizontal and vertical aligned bounding boxes, each row in the
%         format [x y w h]
%
% T is cell array, and each cell contains bounding boxes of the track
% in the format [frame_id x y w h]

[rows cols nfms] = size(fmap);
N = rows * cols;
er_grid_r = min(cols / 50, rows / 50);
% quantize the colors
%tmp = randsample(nfms, 5);
tmp = 1:5;
h = fspecial('disk',5);
for i = 1:size(imgs, 4)
    imgs(:,:,:,i) = imfilter(imgs(:,:,:,i), h);
end

for i = 1:length(tmp)
    cim((i-1)*rows + 1 : i*rows, :, 1) = imgs(:,:,1,tmp(i));
    cim((i-1)*rows + 1 : i*rows, :, 2) = imgs(:,:,2,tmp(i));
    cim((i-1)*rows + 1 : i*rows, :, 3) = imgs(:,:,3,tmp(i));
end

% pre-compute motion and color histogram index maps
[~, cmap] = rgb2ind(lab2uint8(rgb2lab(cim)), params.er_cnbins, 'nodither');
colorMaps = zeros(rows, cols, nfms);
colorPriors = zeros(nfms, params.er_cnbins);
u = zeros(rows, cols, nfms);
v = zeros(rows, cols, nfms);
bu = zeros(rows, cols, nfms);
bv = zeros(rows, cols, nfms);
for i = 1:nfms
    colorMaps(:,:,i) = medfilt2(rgb2ind(lab2uint8(rgb2lab(imgs(:,:,:,i))), cmap) + 1, [5 5], 'symmetric');
    colorPriors(i, :) = histc(colorMaps((i-1)*N+1 : i*N)', 1:params.er_cnbins);
    colorPriors(i, :) = colorPriors(i, :) / sum(colorPriors) + eps;
    u(:,:,i) = medfilt2(flow.u(:,:,i), [3 3], 'symmetric');
    v(:,:,i) = medfilt2(flow.v(:,:,i), [3 3], 'symmetric');
    bu(:,:,i) = medfilt2(flow.bu(:,:,i), [3 3], 'symmetric');
    bv(:,:,i) = medfilt2(flow.bv(:,:,i), [3 3], 'symmetric');
end

% compute tracks
[rootSegs subtrees]  = treeProp2(trees, fmap, imgs, colorMaps, colorPriors, u, v, bu, bv, params);
[C tmap] = treeTrack(fmap, imgs, rootSegs, params);
trackIDS = unique([C(:).tid]);
tmap2 = containers.Map(trackIDS, 1:length(trackIDS));
fids = [C(:).fid];
[~, idx] = sort(fids);
C = C(idx);
T = cell(1, length(trackIDS));
for i = 1:length(C)
    id = tmap2(C(i).tid);
    [ys xs] = ind2sub([rows cols], C(i).PixelIdxList);
    T{id} = [T{id}; C(i).fid min(xs) min(ys) max(xs)-min(xs)+1 max(ys)-min(ys)+1];
end

% compute actionlets
tree_id = [subtrees(:).tree_id];
A = [];
for i = 1:length(trees)
    if params.verbose > 0
        fprintf('Processing frame %d \n', i);
    end
    
    flg = (tree_id == i);
    if ~any(flg)
        continue;
    end
    color_map = colorMaps(:,:,i);
    tree = trees{i};
    fmin = max(1, i - params.er_winsz);
    fmax = min(nfms, i + params.er_winsz);
    dt = treeDecendt(tree);
    st = subtrees(flg);
    for j = 1:length(st)
        if params.verbose > 0
            fprintf('\t subtree %d... \n', j);
        end
        tid = tmap(tmap(:,1) == st(j).track_id, 2);
        if ~isKey(tmap2, tid)
            continue;
        end
        a.tid = tmap2(tid);
        % add the root
        a.start = min(st(j).bbox(:,1));
        a.end = max(st(j).bbox(:,1));
        a.bbox = st(j).bbox(:,2:end);
        a.isRoot = 1;
        a.pid = length(A)+1;
        A = [A a];
        segs = dt{st(j).node_id};
        a.isRoot = 0;
        % propagage every segment of subtree
        for k = 1:length(segs)
            if params.verbose > 0
                fprintf('\t\t part %d: ', k);
            end
            % color prior (keep the same during tracking)
            id = segs(k);
            prm.cp = zeros(1, params.er_cnbins);
            color_counts = histc(color_map(tree(id).PixelIdxList), 1:params.er_cnbins);
            prm.cp = color_counts / sum(color_counts);
            [ys xs] = ind2sub([rows, cols], tree(id).PixelIdxList);
            bbox = [i min(xs) min(ys) max(xs)-min(xs)+1 max(ys)-min(ys)+1];
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
                prm.PixelIdxList = r;
                [ys xs] = ind2sub([rows, cols], prm.PixelIdxList);
                bbox = [bbox; fid min(xs) min(ys) max(xs) - min(xs) + 1  max(ys) - min(ys) + 1];
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
                prm.PixelIdxList = r;
                [ys xs] = ind2sub([rows, cols], prm.PixelIdxList);
                bbox = [fid min(xs) min(ys) max(xs) - min(xs) + 1  max(ys) - min(ys) + 1; bbox];
                if params.verbose > 0
                    fprintf('%d ', fid);
                end
            end
            a.start = min(bbox(:,1));
            a.end = max(bbox(:,1));
            a.bbox = bbox(:,2:end);
            A = [A a];
            if params.verbose > 0
                fprintf('\n');
            end
        end
    end
end
    function seg = trackSegNested
        % TRACKSEG   Using color, motion and spatial priors to track a segment

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


