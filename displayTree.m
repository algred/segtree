function displayTree(img, tree)
% DISPLAYTREE   Draws the bounding boxes of segments in a segment tree

imshow(img);
hold on;
for i = 2:length(tree)
    c = tree(i).BoundingBox;
    plot(c(1,[1:end 1]),c(2,[1:end 1]),'r');
    %text(c(1,1), c(2,1), num2str(i));
end
hold off;