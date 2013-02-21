addpath(genpath(pathstring('Y:\projects\actionlet\v2\external')));
dataroot = pathstring('/research/wvaction/projects/actionlet/v2/demo/1_31');
dataroot2 = pathstring('/research/wvaction/projects/actionlet/v2/demo/2_6');
load(pathstring('Y:\projects\actionlet\matlab\UCFSports\ucf_sports_data.mat'));
outroot = pathstring('/research/wvaction/projects/actionlet/v2/demo/2_14');



for cls = 8:length(class_names)
    idx = find(class_labels == cls);
    for i = 3:3
        vid = idx(i);
        % %         wvname = [outroot filesep  class_names{class_labels(vid)} '_' video_list(vid).vname '_root.avi'];
        % %         if exist(wvname, 'file')
        % %             continue;
        % %         end
        % %         wobj = VideoWriter(wvname, 'Uncompressed AVI');
        % %         wobj.FrameRate = 1;
        % %         open(wobj);
        
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
        %         %[C S] = propCandSegs(tree, fmap, imgs, flow, M, estRootParams);
                 load([outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_cand.mat']);
        CC = getRootCandidate(C, fmap, imgs, estRootParams);
        [R S] = estimateRoot(CC, estRootParams);
        fid = [CC(:).fid];
        cmap = colormap(lines(length(R)));
        colormap(jet);
        for k = 1:size(M, 3)
            C1 = CC(fid == k);
            X = zeros(length(C1) * 2, 4);
            img1 = imgs(:,:,:,k);
            img2 = zeros(rows, cols, 3);
            img3 = img2;
            bgidx = 1:N;
            for j = 1:length(C1)
                [ys xs] = ind2sub([rows, cols], C1(j).PixelIdxList);
                if size(xs, 1) > 1
                    xs = xs';
                    ys = ys';
                end
                X((j-1)*2 + 1: j*2, :) = boundingBox(xs, ys);
                bgidx = setdiff(bgidx, C1(j).PixelIdxList);
            end
            img3 = imgs(:,:,:,k);
            [~, bestTrackId] = max(S);
            T = R(bestTrackId);
            for kk = 1:length(T)
                rfid = [T{kk}.fid];
                R1 = T{kk}(rfid == k);
                if isempty(R1)
                    continue;
                end
                if length(R1) >= 2
                    fprintf('ERROR: More than one matching box on a single frame\n');
                    keyboard;
                end
                img2(R1.PixelIdxList) = img1(R1.PixelIdxList)/2 + cmap(kk, 1)/2;
                img2(R1.PixelIdxList + N) = img1(R1.PixelIdxList + N)/2 + cmap(kk, 2)/2;
                img2(R1.PixelIdxList + 2*N) = img1(R1.PixelIdxList + 2*N)/2 + cmap(kk, 3)/2;
                img2 = uint8(min(255, max(1, round(img2))));
                [ys xs] = ind2sub([rows, cols], R1.PixelIdxList);
                if size(xs, 1) > 1
                    xs = xs';
                    ys = ys';
                end
                img3 = drawBox(img3, boundingBox(xs, ys), ceil(254 * cmap(kk, :)));
            end
            
            
            img1(bgidx) = 0;
            img1(bgidx + N) = 0;
            img1(bgidx + 2 * N) = 0;
            
            subplot(3, 2, 2); imagesc(fmap(:,:,k));
            subplot(3, 2, 1); imagesc(drawTree(imgs(:,:,:,k), tree{k}));
            subplot(3, 2, 3); imagesc(drawBox(imgs(:,:,:,k), X, [255, 0, 0]));
            subplot(3, 2, 4); imagesc(img1);
            subplot(3, 2, 5); imagesc(img3);
            subplot(3, 2, 6); imagesc(img2);
            keyboard;
            % %             F = getframe(gcf);
            % %             writeVideo(wobj, F.cdata);
        end
        fprintf('Finished video %s\n', video_list(vid).vname);
        % %         save([outroot filesep class_names{class_labels(vid)} '_' video_list(vid).vname '_cand.mat'], 'C', 'S', 'R');
        % %         close(wobj);
    end
end
