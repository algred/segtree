function seg = trackSeg(prm)
% TRACKSEG   Using color, motion and spatial priors to track a segment
%
% fields of prm:
%   cp, mp, color_map, motion_map, x0, y0, sp_sz, PixelIdxList
%

% compute spatial prior
try
    
[rows cols] = size(prm.color_map);
h = fspecial('gaussian', prm.sp_sz * 4, prm.sp_sz);
spatial_pm = zeros(rows, cols);
radius = round(prm.sp_sz/2);
xidx = max(1, prm.x0 - radius) : min(cols, prm.x0 + radius);
yidx = max(1, prm.y0 - radius) : min(rows, prm.y0 + radius);
spatial_pm(yidx,  xidx) = h(yidx - prm.y0 + prm.sp_sz * 2, xidx - prm.x0 + prm.sp_sz * 2);
% motion prior
motion_pm = ni(medfilt2(reshape(prm.mp(prm.motion_map(:)), rows, cols), [5 5], 'symmetric'));
% color prior
color_pm = zeros(rows, cols);
for i = 1:size(prm.cp, 1)
    cp = prm.cp(i,:);
    color_pm = color_pm + medfilt2(reshape(cp(prm.color_map(:)), rows, cols), [5 5], 'symmetric');
end
color_pm = ni(color_pm);

% find the match
pm = ni(color_pm.*motion_pm.*spatial_pm);
pm2 = round(pm * 20);
tsz = length(prm.PixelIdxList);
minszdiff = inf;
for i = 0:20
    bw = imclose((pm2 >= i), prm.se);
    CC = bwconncomp(bw);
    for j = 1:CC.NumObjects
        szdiff = abs(length(CC.PixelIdxList{j}) - tsz);
        if szdiff < minszdiff
            minszdiff = szdiff;
            seg = CC.PixelIdxList{j};
        end
    end
end

tsz1 = length(seg);
if tsz1 <= tsz * 0.3 || tsz1 >= tsz * 1.3
    seg = [];
    return;
end

catch exception
    getReport(exception)
    keyboard;
end

end