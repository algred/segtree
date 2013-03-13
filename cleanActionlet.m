function A = cleanActionlet(A, min_len)
% CLEANACTIONLET    Clean actionlets in actionlet array A

% remove short actionlets
S = [A(:).start];
E = [A(:).end];


D = ((E - S + 1) > min_len);
keep = find(D);
removed = setdiff(1:length(A), keep);

P = zeros(size(A));
P(keep) = 1:length(keep);
P(removed) = ones(size(removed))*-1;

A = A(D);
P2 = [A(:).pid]';
P2 = num2cell(P(P2), 1);
[A(:).pid] = P2{:};

end