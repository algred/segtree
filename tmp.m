%load(pathstring('Y:\projects\actionlet\matlab\HOHA1\hoha1_data.mat'));

load(pathstring('Y:\projects\actionlet\v2\data\ucf_sports\ucf_sports_data.mat'));

addpath(genpath(pathstring('Y:\projects\actionlet\v2\external')));
addpath(genpath(pathstring('Y:\projects\actionlet\v2\segtree')));
addpath(genpath(pathstring('Y:\projects\actionlet\v2\hoha1')));
addpath(genpath(pathstring('Y:\projects\actionlet\v2\ucfsports')));

flowpath = pathstring('Y:\data\video_data\features\brox_flow\hoha1\segments');
treepath = pathstring('Y:\projects\actionlet\v2\data\ucf_sports\segtree');
if ispc
    imgs = read(pathstring(VideoReader(video_list(vid).videoname)));
else
    imgs = read_video(pathstring(video_list(vid).videoname), 1);
end

% load([flowpath filesep num2str(vid) '_flow.mat']);
% load([treepath filesep num2str(vid) '_trees.mat']);
% 
% nfms = min(size(imgs, 4), video_list(vid).end);
% imgs = imgs(:,:,:, video_list(vid).start: nfms);

load(pathstring(video_list(vid).flowname));
load([treepath filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_trees.mat']);

rows = size(imgs, 1);
cols = size(imgs, 2);

params = er_params_hoha(rows, cols);
params.alpha_AB = 1.9;
params.alpha_M = 2;


for frame_id = 1 : size(imgs, 4)-1
    I = imgs(:,:,:,frame_id);
    subplot(3, 3, 1); imshow(I);
    subplot(3, 3, 2); imagesc(gb_CSMGS(:,:,frame_id)); colorbar;
    subplot(3, 3, 3); imagesc(ucms(:,:,frame_id)); colorbar;
    subplot(3, 3, 4); imagesc(fmap(:,:,frame_id)); colorbar;
    subplot(3, 3, 5); displayTree(I, trees3{frame_id});
    subplot(3, 3, 6); displayTree(I, trees5{frame_id});
    u = flow.u(:,:,frame_id); v = flow.v(:,:,frame_id); 
    hu = flow.hu(:,:,frame_id); hv = flow.hv(:,:,frame_id);
    [tree tree2 tree3 tree4 ucm gb_CSMG] = segtree(imgs(:,:,:,frame_id), u, v, hu, hv, params);
    subplot(3, 3, 7); imagesc(gb_CSMG);
    subplot(3, 3, 8); imagesc(ucm);
    subplot(3, 3, 9); displayTree(I, tree4);
    keyboard;
end