addpath(genpath(pathstring('Y:\projects\actionlet\v2\external')));
dataroot = pathstring('/research/wvaction/projects/actionlet/v2/demo/1_31');
load(pathstring('Y:\projects\actionlet\matlab\UCFSports\ucf_sports_data.mat'));
outroot = pathstring('/research/wvaction/projects/actionlet/v2/demo/2_6');

params.pr_maxbgm = 30;
params.pr_cnbins = 128;
params.minsz_ratio = 0.003;
params.min_motion = 0.8;
params.maxsz_ratio = 0.25;
params.min_shrink_ratio = 0.67;
params.DISPLAY = 0;
params.fg_ratio = 0.3;
params.bw_quant_num = 5;
params.bw_sigma_blur = 30;

for cls = 1:length(class_names)
    idx = find(class_labels == cls);
    for i = 3:3
        vid = idx(i);        
        load([dataroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_trees.mat']);
        load(pathstring(video_list(vid).flowname));
        M = motionMagnitude(flow);
        if ispc 
            imgs = read(VideoReader(pathstring(video_list(vid).videoname)));
        else
            imgs = read_video(pathstring(video_list(vid).videoname), 1);
        end
        [cp mp cpmap mpmap fmap] = estimatePriors(tree, imgs, M, params);
        for fid = 1:size(fmap, 3)
            [t1 t2] = segtree_with_prior(imgs(:,:,:,fid), flow.u(:,:,fid), flow.v(:,:,fid), flow.hu(:,:,fid), flow.hv(:,:,fid), fmap(:,:,fid), params);
            org_tree2{fid} = t1;
            tree2{fid} = t2;
        end
        
        colormap(jet);
        for k = 1:size(M, 3)
            subplot(2, 2, 1); imshow(imgs(:,:,:,k));
            subplot(2, 2, 2); imshow(drawTree(imgs(:,:,:,k), tree{k}));
            subplot(2, 2, 3); imagesc(fmap(:,:,k)); colorbar;
            subplot(2, 2, 4); imshow(drawTree(imgs(:,:,:,k), tree2{k}));
            keyboard;
        end
    end
end
