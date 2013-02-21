outroot = '/research/wvaction/projects/actionlet/v2/demo/2_19';

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


for cls = 1:length(class_names)
    idx = find(class_labels == cls);
    for i = 1:3
        vid = idx(i);
        wvname = [outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_segment.avi'];
        if exist('wvname', 'file')
            continue;
        end
        clear('tree', 'org_tree');
        imgs = read_video(pathstring(video_list(vid).videoname), 1);
        load(pathstring(video_list(vid).flowname));
        rows = size(imgs, 1);
        cols = size(imgs, 2);
           
        wobj = VideoWriter(wvname, 'Uncompressed AVI');
        wobj.FrameRate = 1;
        open(wobj);
        
        for frame_id = 1 : size(imgs, 4)-1
            I = imgs(:,:,:,frame_id);
            [tree1 tree2] = segtree(I, flow.u(:,:,frame_id),flow.v(:,:,frame_id),...
                flow.hu(:,:,frame_id), flow.hv(:,:,frame_id), params);
            tree{frame_id} = tree2;
            org_tree{frame_id} = tree1;
            img1 = drawTree(imgs(:,:,:,frame_id), tree2);
            %imshow(img1);
            writeVideo(wobj, img1);
        end
        
        save([outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_trees.mat'], 'tree', 'org_tree');
        close(wobj);
    end
end