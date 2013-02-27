M = zeros(nt);
for i = 1:nt-1
    mint = min(D{i});
    maxt = max(D{i});
    for j = i+1:nt
        if mint > max(D{j})
            M(i, j) = getSimilarity(C(I(i,1)).Feat, C(I(j,2)).Feat, 1:3, params.er_cov_disturb, params.er_cov_scale);
        elseif maxt < min(D{j})
            M(i, j) = getSimilarity(C(I(i,2)).Feat, C(I(j,1)).Feat, 1:3, params.er_cov_disturb, params.er_cov_scale);
        end
        M(i, j) = M(j, i);
    end
end