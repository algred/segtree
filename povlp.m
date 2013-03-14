function ratio = povlp(r1, r2)
% POVLP   computes the overlap between the corresponding rectangle in r1
% and r2. r1 and r2 should have the same number of rows
% 
% r1, r2: input rectangles in the form [x y w h]
% ratio: ratios between overlapped area and their respective area

n = size(r2, 1);
x1 = max(r1(:, 1), r2(:,1));
y1 = max(r1(:, 2), r2(:,2));
x2 = min(r1(:, 1) + r1(:, 3)-1, r2(:,1) + r2(:,3)-1);
y2 = min(r1(:, 2) + r1(:, 4)-1, r2(:,2) + r2(:,4)-1);
a1 = r1(:, 3).*r1(:, 4);
a2 = r2(:, 3).*r2(:, 4);
ratio = zeros(n, 2);
flg = ((x2 - x1) > 0 & (y2 - y1) > 0);
if ~any(flg)
    return;
end
ovlp = (x2(flg)-x1(flg)+1).*(y2(flg)-y1(flg)+1);
ratio(flg, 1) = ovlp ./ a1(flg);
ratio(flg, 2) = ovlp ./ a2(flg);