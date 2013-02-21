addpath(genpath(pathstring('Y:\projects\actionlet\v2\external')));
dataroot = pathstring('/research/wvaction/projects/actionlet/v2/demo/1_31');
dataroot2 = pathstring('/research/wvaction/projects/actionlet/v2/demo/2_6');
load(pathstring('Y:\projects\actionlet\matlab\UCFSports\ucf_sports_data.mat'));
outroot = pathstring('/research/wvaction/projects/actionlet/v2/demo/2_14');


 
for cls = 3:length(class_names)
    idx = find(class_labels == cls);
    for i = 3:3
        vid = idx(i);
        load([dataroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_trees.mat']);
        load([dataroot2 filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_priors.mat']);
        load(pathstring(video_list(vid).flowname));
        M = motionMagnitude(flow);
        if ispc
            imgs = read(VideoReader(pathstring(video_list(vid).videoname)));
        else
            imgs = read_video(pathstring(video_list(vid).videoname), 1);
        end
        rows = size(imgs, 1);
        cols = size(imgs, 2);
        N = rows * cols;
        estRootParams = er_params(rows, cols);
        estRootParams.verbose = 1;
        tic;
        [rootSegs trackID] = treeProp(tree, fmap, imgs, flow, M, estRootParams);
        toc;
        tracks = unique(trackID);
        cmap = colormap(lines(length(tracks)));
        colormap(jet);
        
        for k = 1:size(M, 3)
            thisRootSeg = rootSegs{k};
            im1 = imgs(:,:,:,k);
            im2 = im1;
            fidx = [];
            for h = 1:length(thisRootSeg)
                color_id = find(tracks == thisRootSeg(h).TemplateId);
                [ys xs] = ind2sub([rows cols], thisRootSeg(h).PixelIdxList);
                fidx = union(fidx, thisRootSeg(h).PixelIdxList);
                im2 = drawBox(im2, boundingBox(xs, ys), ceil(254 * cmap(color_id, :)));
            end
            bgidx = setdiff(1:N, fidx);
            im1(bgidx) = 0;
            im1(bgidx + N) = 0;
            im1(bgidx + 2 * N) = 0;
            
            subplot(2,1,1); imshow(im1);
            subplot(2,1,2); imshow(im2);
            keyboard;
        end
        
    end
    fprintf('Finished video %s\n', video_list(vid).vname);
end

