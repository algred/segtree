function [g L] = trackAnnotate(A, T, bboxs, params)

g = zeros(size(A));
for i = 1:length(A)
    fid = A(i).start:A(i).end;
    boxs1 = A(i).bbox;
    for j = 1:length(bboxs)
        boxs2 = bboxs{j}(fid, :);
        flg = (boxs2(:, 1) > 0);
        if ~any(flg)
            break;
        end
        ovlp = povlp(boxs1(flg, :), boxs2(flg, :));
        if A(i).isRoot > 0
            idx = (ovlp(:,1) > params.train_ovlp) & (ovlp(:,2) > params.train_ovlp);
            if any(idx)
                g(i) = j;
                break;
            end  
        else
            idx = ovlp(:,1) > params.train_ovlp;
            if any(idx)
                g(i) = j;
                break;
            end  
        end
    end
end


D = [A(:).tid];
R = ([A(:).isRoot] > 0);
P = ~R;
L = zeros(size(T));
for i = 1:length(T)
    flg = (D == i);
    if any(g(flg & R) > 0) && any(g(flg & P) > 0)
        L(i) = 1;
    else
        g(flg) = -1;
    end
end