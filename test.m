video_id = 15;
frame_off = 1;
% % % 
addpath(genpath(pathstring('Y:\projects\actionlet\v2\external')));
load(pathstring('Y:\projects\actionlet\matlab\UCFSports\ucf_sports_data.mat'));
% %  
imgs = read_video(pathstring(video_list(video_id).videoname), 1);
load(pathstring(video_list(video_id).flowname));
rows = size(imgs, 1);
cols = size(imgs, 2);

% % wobj = VideoWriter([class_names{class_labels(video_id)} '_' video_list(video_id).vname '_segment.avi'], 'Uncompressed AVI');
% % wobj.FrameRate = 1;
% % open(wobj);
% 
 params.minsz_ratio = 0.003;
 params.min_motion = 0.8;
 params.maxsz_ratio = 0.25;
 params.min_shrink_ratio = 0.67;
 params.DISPLAY = 1;
 params.fg_ratio = 0.3;
 params.bw_quant_num = 5;
 params.bw_sigma_blur = 30;
 params.cov_disturb = 1e-6;
 params.cov_scale = 2;
 
 for frame_id = 1: size(imgs, 4)-1
     I = imgs(:,:,:,frame_id);
     [tree1 tree2] = segtree(I, flow.u(:,:,frame_id),flow.v(:,:,frame_id),...
          flow.hu(:,:,frame_id), flow.hv(:,:,frame_id), params); 
     tree{frame_id} = tree2;
     org_tree{frame_id} = tree1;
     keyboard;
%     img1 = drawTree(imgs(:,:,:,frame_id), tree2);
%     imshow(img1);    
%     h = figure('Name', 'Tree'); 
%     set(0,'CurrentFigure',h)
%     displayTree2(I, tree{frame_id}, 1);
 end

% nfeat = size(tree{frame_off + 1}(2).Feat, 1);
% for i = frame_off + 1:frame_off + 4%length(tree)-1
%    S1{i} = treeSim(tree{i}, tree{i+1}, [6:10], params);
% end

