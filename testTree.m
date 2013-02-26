outroot = pathstring('/research/wvaction/projects/actionlet/v2/demo/2_25');

addpath(genpath(pathstring('Y:\projects\actionlet\v2\external')));
addpath(pathstring('Y:\code\bsr\lib'));
load(pathstring('Y:\projects\actionlet\matlab\UCFSports\ucf_sports_data.mat'));

params.minsz_ratio = 0.003;
params.min_motion = 0.8;
params.maxsz_ratio = 0.25;
params.min_shrink_ratio = 0.67;
params.DISPLAY = 0;
params.fg_ratio = 0.3;
params.bw_quant_num = 5;
params.bw_sigma_blur = 30;
params.pr_maxbgm = 30;
params.pr_cnbins = 128;
params.er_score_thre = 0.2;
params.er_delta = 3;

for cls = 2:length(class_names)
    idx = find(class_labels == cls);
    for i = 3:3
        vid = idx(i);
        wvname = [outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_tree.avi'];
        if exist('wvname', 'file')
            continue;
        end
        wobj = VideoWriter(wvname, 'Uncompressed AVI');
        wobj.FrameRate = 1;
        open(wobj);
        
        clear('trees', 'trees2', 'trees3', 'trees4', 'ucms', 'gb_CSMGS');
          imgs = read_video(pathstring(video_list(vid).videoname), 1);
        load(pathstring(video_list(vid).flowname));
        M = motionMagnitude(flow);
        rows = size(imgs, 1);
        cols = size(imgs, 2);

        for frame_id = 1 : size(imgs, 4)-1
            I = imgs(:,:,:,frame_id);
            [tree tree2 tree3 tree4 ucm gb_CSMG] = segtree(I, flow.u(:,:,frame_id),flow.v(:,:,frame_id),...
                flow.hu(:,:,frame_id), flow.hv(:,:,frame_id), params);
            trees{frame_id} = tree;
            trees2{frame_id} = tree2;
            trees3{frame_id} = tree3;
            trees4{frame_id} = tree4;
            ucms{frame_id} = ucm;
            gb_CSMGS{frame_id} = gb_CSMG;
        end
        [cp cpmap mpmap fmap] = estimatePriors2(trees4, imgs, M, params);
        trees5 = rmBackground(trees4, fmap, params);
        for frame_id = 1 : size(imgs, 4)-1
            I = imgs(:,:,:,frame_id);
            subplot(3, 3, 1); imshow(I);
            subplot(3, 3, 2); imagesc(gb_CSMGS{frame_id}); colorbar;
            subplot(3, 3, 3); imagesc(ucms{frame_id}); colorbar;
            subplot(3, 3, 4); displayTree(I, trees{frame_id});
            subplot(3, 3, 5); displayTree(I, trees2{frame_id});
            subplot(3, 3, 6); displayTree(I, trees3{frame_id});
            subplot(3, 3, 7); displayTree(I, trees4{frame_id});
            subplot(3, 3, 8); imagesc(fmap(:,:,frame_id));
            subplot(3, 3, 9); displayTree(I, trees5{frame_id});
            keyboard;
            frame = getframe(gcf);
            writeVideo(wobj, frame);
        end
        
        save([outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_trees.mat'], 'trees5', 'trees3', 'ucms','gb_CSMGS', 'fmap');
        close(wobj);
    end
end