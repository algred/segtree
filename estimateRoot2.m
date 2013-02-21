function [R S] = estimateRoot(C, params)
% ESTIMATEROOT    Estimate the human area given the foreground 
% estimation and candidate positions
%
% Input
% fmap: foreground probability maps
% imgs: video frames
% C: a structure array containing candidate regions 
%
% Output:
% R: a cell array. Each cell contains a structure array which may
% correspond to one track of a human 


% initialize the tracks
fid = [C(:).fid];
nfms = max(fid);
cidx = find(fid == 1);
for i = 1:length(cidx)
    I{i} = cidx(i);
    S(i) = C(cidx(i)).Score;
end

% grow the tracks
for i = 2:nfms
    cidx = find(fid == i);
    [~,tidx] = sort(S, 'descend');
    ntracks = length(tidx);
    for j = 1:ntracks
        seg = C(I{tidx(j)}(end));
        sim = zeros(size(cidx));
        for k = 1:length(cidx)
            sim(k) = getSimilarity(seg.Feat, C(cidx(k)).Feat, ...
                1:3, params.er_cov_disturb, params.er_cov_scale);
        end
        [maxsim id] = max(sim);
        if maxsim > params.er_sim_thre
           I{tidx(j)} = [I{tidx(j)} cidx(id)];
           S(tidx(j)) = S(tidx(j)) + C(cidx(id)).Score;
           cidx = cidx([1:id-1 id+1:end]);
           if isempty(cidx)
               break;
           end
        end
    end
    for j = 1:length(cidx)
        I{ntracks + 1} = cidx(j);
        S(ntracks + 1) = C(cidx(j)).Score;
    end
end

% construct the result
for i = 1:length(I)
    R{i} = C(I{i});
end

end