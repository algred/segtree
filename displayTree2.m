function displayTree2(img, tree, DISPLAY_BB)
% DISPLAYTREE   Draws the bounding boxes of segments in a segment tree

[rows cols nch] = size(img);
N = rows * cols;
tree = treeDepth(tree);

D = [tree(:).Depth];
P = [tree(:).Parent];
childs = treeDecendt(tree);

if ~exist('DISPLAY_BB', 'var')
    DISPLAY_BB = 0;
end

for i = 2:length(tree)
    try
    % special null node
    if isempty(tree(i).PixelIdxList)
        segments{i} = uint8(zeros(30, 30, 3));
        continue;
    end
    
    % normal nodes
    c = tree(i).BoundingBox;
    c(1,:) = min(cols, max(1, round(c(1,:))));
    c(2,:) = min(rows, max(1, round(c(2,:))));
    minx = min(c(1,:));
    miny = min(c(2,:));
    maxx = max(c(1,:));
    maxy = max(c(2,:));
    if DISPLAY_BB
        im = uint8(zeros(rows, cols, 3));
        bw = poly2mask(c(1, [1:end 1]), c(2, [1:end 1]), rows, cols);
        idxList = find(bw);
        idxList2 = [idxList; idxList+N; idxList+2*N];
        im(idxList2) = img(idxList2);
    else
        im = uint8(ones(rows, cols, 3) * 4);
        im(tree(i).PixelIdxList) = img(tree(i).PixelIdxList);
        im(tree(i).PixelIdxList + N) = img(tree(i).PixelIdxList + N);
        im(tree(i).PixelIdxList + 2*N) = img(tree(i).PixelIdxList + 2*N);
    end
    
    % draw the orientation line
    [ys xs] = ind2sub([rows cols], tree(i).PixelIdxList);
    y1 = min(ys);
    y2 = max(ys);
    if tree(i).or == pi/2
        ys = [y1:y2]';
        xs = ones(size(ys))*tree(i).x0;
    elseif tree(i).or == 0
        xs = [min(xs):max(xs)]';
        ys = ones(size(xs))*tree(i).y0;
    else
        x1 = tree(i).x0 + abs(y1 - tree(i).y0) / tan(tree(i).or);
        x2 = tree(i).x0 - abs(y2 - tree(i).y0) / tan(tree(i).or);
        [xs ys] = bresenham(x1,y1,x2,y2);
        flg = (xs >= 1 & xs <= cols);
        xs = xs(flg);
        ys = ys(flg);
    end
    p = sub2ind([rows cols], ys, xs);
    im(p) = 255;
    im(p+N) = 0;
    im(p+2*N) = 0;
    segments{i} = imcrop(im, [minx miny maxx-minx+1 maxy-miny+1]);
    catch exception
        getReport(exception)
        keyboard;
    end
end



plotrow = max(D);

C = zeros(length(tree), plotrow);
for i = length(tree):-1:1
    child = childs{i};
    if isempty(child)
        for j = plotrow:-1:D(i)
            C(i, j) = 1;
        end
        continue;
    end
    for j = plotrow:-1:D(i)
        C(i, j) = sum(C(P == i, j));
    end
end



% for i = 1:length(tree)
%     child = childs{i};
%     if isempty(child)
%         continue;
%     end
%     for j = 1:plotrow
%         C(i, j) = sum(D(child) == j);
%     end
% end

maxW = max(1, max(C, [], 2));
plotcol = 0;
for i = 1:max(D)
    val = sum(maxW(D == i));
    if plotcol < val
        plotcol = val;
    end
end

offset = zeros(plotrow, 1);

% depth first
D = D - 1;
stack = find(D == 1);
while ~isempty(stack)
    id = stack(end);
    stack = stack(1:end-1);
    depth = D(id);
    child = childs{id};
    subplot(plotrow, plotcol, (depth - 1) * plotcol + offset(depth) + 1); imshow(segments{id});
    offset(depth + 1) = offset(depth);
    offset(depth) = offset(depth) + maxW(id);
    if isempty(child)
        continue;
    end
    stack = [stack child( D(child) == depth + 1 )];
end






