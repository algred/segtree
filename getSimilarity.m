% return the simliarty score between two tree nodes
function s=getSimilarity(a,b,chs, disturb, scale)

try
cov_a=cov(a(chs,:)');
cov_b=cov(b(chs,:)');

%for numeric stability

[U_a,S_a,V_a]=svd(cov_a);
v = diag(S_a);
v = v + disturb./(v+eps);
S_a2 = diag(v);
cov_a = U_a * S_a2 * V_a';

[U_b,S_b,V_b]=svd(cov_b);
v = diag(S_b);
v = v + disturb./(v+eps);
S_b2 = diag(v);
cov_b = U_b * S_b2 * V_b';

e=abs(eig(cov_a,cov_b)+eps);

s=sum(log(e).^2);

s = exp(-s/scale);

if isnan(s)
    keyboard;
end

catch exception
    getReport(exception)
    keyboard;
end