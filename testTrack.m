outroot = pathstring('/research/wvaction/projects/actionlet/v2/demo/2_25');

addpath(genpath(pathstring('Y:\projects\actionlet\v2\external')));
addpath(pathstring('Y:\code\bsr\lib'));
load(pathstring('Y:\projects\actionlet\matlab\UCFSports\ucf_sports_data.mat'));


for cls = 8:length(class_names)
    idx = find(class_labels == cls);
    for i = 3:3
        vid = idx(i);
        wvname = [outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_track.avi'];
        if exist(wvname, 'file')
            continue;
        end
        wobj = VideoWriter(wvname, 'Uncompressed AVI');
        wobj.FrameRate = 1;
        open(wobj);
        
        imgs = read_video(pathstring(video_list(vid).videoname), 1);
        load(pathstring(video_list(vid).flowname));
        M = motionMagnitude(flow);
        rows = size(imgs, 1);
        cols = size(imgs, 2);
        params = er_params(rows, cols);
        
        load([outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_trees.mat']);
        [rootSegs subtrees]  = treeProp(trees5, fmap, imgs, flow, M, params);
          [C tmap] = treeTrack(fmap, imgs, rootSegs, params);
        
        F = [C(:).fid];
        T = [C(:).tid];
        tids = unique(T);
        cmap = colormap(lines(length(tids)));
        for frame_id = 1 : size(imgs, 4)-1
            I = imgs(:,:,:,frame_id);
            subplot(2, 1, 1); displayTree(I, trees5{frame_id});
            C1 = C(F == frame_id);
            subplot(2, 1, 2); imshow(I);
            hold on;
            for h = 1:length(C1)
                color_id = find(tids == C1(h).tid);
                [ys xs] = ind2sub([rows cols], C1(h).PixelIdxList);
                c = boundingBox(xs, ys);
                plot(c(1,[1:end 1])',c(2,[1:end 1])','Color', cmap(color_id, :));
                text(c(1,1), c(2,1), num2str(C1(h).tid), 'FontSize', 20);
            end
            hold off;
            keyboard;
            frame = getframe(gcf);
            writeVideo(wobj, frame);
        end
        save([outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_tracks.mat'], 'C', 'subtrees', 'rootSegs','tmap');
        close(wobj);
        fprintf('over\n');
        keyboard;
    end
end
