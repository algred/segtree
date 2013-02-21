addpath(genpath(pathstring('Y:\projects\actionlet\v2\external')));
dataroot = pathstring('/research/wvaction/projects/actionlet/v2/demo/1_31');
load(pathstring('Y:\projects\actionlet\matlab\UCFSports\ucf_sports_data.mat'));
outroot = pathstring('/research/wvaction/projects/actionlet/v2/demo/2_6');

params.pr_maxbgm = 30;
params.pr_cnbins = 128;

for cls = 1:length(class_names)
    idx = find(class_labels == cls);
    for i = 3:3
        vid = idx(i);
%         wvname = [outroot filesep  class_names{class_labels(vid)} '_' video_list(vid).vname '_priors.avi'];
%         if exist(wvname, 'file')
%             continue;
%         end
%         wobj = VideoWriter(wvname, 'Uncompressed AVI');
%         wobj.FrameRate = 1;
%         open(wobj);
        
        load([dataroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_trees.mat']);
        old_fmap = load([outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_priors.mat']);
        load(pathstring(video_list(vid).flowname));
        M = motionMagnitude(flow);
        if ispc 
            imgs = read(VideoReader(pathstring(video_list(vid).videoname)));
        else
            imgs = read_video(pathstring(video_list(vid).videoname), 1);
        end
        [cpmap mpmap fmap] = estimatePriors3(tree, imgs, M, params);
        colormap(jet);
        for k = 1:size(M, 3)
            subplot(3, 2, 1); imshow(imgs(:,:,:,k));
            subplot(3, 2, 2); imshow(drawTree(imgs(:,:,:,k), tree{k}));
            subplot(3, 2, 3); imagesc(cpmap(:,:,k)); colorbar;
            subplot(3, 2, 4); imagesc(mpmap(:,:,k)); colorbar;
            subplot(3, 2, 5); imagesc(old_fmap.fmap(:,:,k)); colorbar;
            subplot(3, 2, 6); imagesc(fmap(:,:,k)); colorbar;
            keyboard;
%             F = getframe(gcf);
%             writeVideo(wobj, F.cdata);
%             fmap(:,:,k) = this_fmap;
        end
%         save([outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_priors.mat'], ...
%             'cp', 'mp', 'cpmap', 'mpmap', 'fmap');
%         close(wobj);
    end
end
