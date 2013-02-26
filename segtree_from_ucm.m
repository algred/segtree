function [tree] = segtree_from_ucm(ucm, params)
% SEGTREE create a hierachical segmentation tree from a UCM

[rows cols] = size(ucm);
ucm = round(ucm * 100);
bryvs = setdiff(unique(ucm(:)), [0:15]);

r.PixelIdxList = 1:cols*rows;
r.x0 = cols / 2;
r.y0 = rows / 2;
r.or = 0;
r.xr = [];
r.yr = [];

r.BoundingBox = [1 cols cols 1; 1 1 rows rows];
r.PixelIdxListBB = 1:cols*rows;
r.Parent = 1;

tree(1) = r;

if length(bryvs) <= 1
    return
end

minsz = params.minsz_ratio * rows * cols;

% We do a breath first traversal
queue = [1];
bryid = length(bryvs);
if params.DISPLAY
    h = figure('Name', 'Segmentation');
    set(0,'CurrentFigure',h);
end
while length(queue) >= 0 && bryid >= 1
    label = bwlabel(ucm <= bryvs(bryid));
    new_queue = [];
    for i = 1:length(queue)
        try
            r = tree(queue(i));
            segmap = label(r.PixelIdxList);
            lv = setdiff(unique(segmap), 0);
            sz1 = length(r.PixelIdxList);
            for j = 1:length(lv)
                pixelList = r.PixelIdxList(segmap == lv(j));
                sz = length(pixelList);
                % not segmented by the current threshold
                if sz1 - sz < minsz
                    new_queue = [new_queue queue(i)];
                    break;
                end
                % ignore small segments
                if sz < minsz
                    continue;
                end
               
                newR.PixelIdxList = pixelList;

                
                % Moments of the region: center and orientation
                [ys xs] = ind2sub([rows cols], newR.PixelIdxList);
                [x0 y0 or] = segProperty([xs' ys']);
                newR.x0 = x0;
                newR.y0 = y0;
                newR.or = or;
                
                % relative coordinates
                cosr = cos(newR.or);
                sinr = sin(newR.or);
                newR.xr = (xs - x0) * cosr + sinr * (ys - y0);
                newR.yr = -(xs - x0) * sinr + cosr * (ys - y0);
                
                % Orientation parallel bounding box
                c = orParallelBB([xs' ys'],  newR.or);
                newR.BoundingBox = c;
                bwBB = poly2mask(c(1,[1:end 1]), c(2, [1:end 1]), rows, cols);
                newR.PixelIdxListBB = find(bwBB);
                newR.Parent = queue(i);
         
                len = length(tree);
                tree(len+1) = newR;
                new_queue = [new_queue len+1];
            end
        catch exception
            getReport(exception)
            keyboard;
        end
    end
    if params.DISPLAY > 0
        fprintf('bryv: %f\n', bryvs(bryid));
        map = zeros(rows, cols);
        for i = 1:length(tree)
            map(tree(i).PixelIdxList) = i;
        end
        %imagesc(map);
        clf;
        imshow(map, lines);
        hold on;
        for i = 2:length(tree)
            c = tree(i).BoundingBox;
            plot(c(1,[1:end 1]),c(2,[1:end 1]),'r');
        end
        hold off;
        keyboard;
    end
    bryid = bryid - 1;
    queue = new_queue;
end


