function A = actionletPos(A, T)

for i = 1:length(A)
    if A(i).isRoot > 0  % just to make the struct array valid
        A(i).x = 0;
        A(i).y = 0;
        continue;
    end
    fid = A(i).start : A(i).end;
    x = zeros(size(fid)); y = zeros(size(fid));
    for j = 1:length(fid)
        bbox1 = A(i).bbox(j, :);
        [~, fid2] = min(abs(T{A(i).tid}(:, 1) - fid(j)));
        bbox2 = T{A(i).tid}(fid2, 2:end);
        x(j) = ((bbox1(1) + bbox1(1) + bbox1(3) - 1) / 2 - (bbox2(1) + bbox2(1) + bbox2(3) - 1) / 2) / bbox2(3);
        y(j) = ((bbox1(2) + bbox1(2) + bbox1(4) - 1) / 2 - (bbox2(2) + bbox2(2) + bbox2(4) - 1) / 2) / bbox2(4);
    end
    A(i).x = x;
    A(i).y = y;
end

end