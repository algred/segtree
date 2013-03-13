function A1 = extendActionlet(A, len_thre, shift)

D = [A(:).end] - [A(:).start] + 1;
idx = find(D > len_thre);

A1 = A;

for i = 1:length(idx)
    id = idx(i);
    a = A(id);
    a1 = a;
    a2 = a;
    a1.end = a.end - shift;
    a2.start = a.start + shift;
    a1.bbox = a1.bbox(1:(end-shift), :);
    a2.bbox = a2.bbox((shift + 1):end, :);
    A1 = [A1 a1 a2];
end

end