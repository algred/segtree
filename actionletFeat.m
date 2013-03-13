function F = actionletFeat(A, imgs, flow, params)
% ACTIONLETFEAT   Compute the features of an array of actionlets 
% F(i,:) contains the feature of actionlet A(i)

[rows cols nfms] = size(flow.u);
hogIntIm = zeros(rows, cols, 9, nfms);
hofIntIm = hogIntIm;
mbhxIntIm = hogIntIm;
mbhyIntIm = hogIntIm;

% Compute the integral images first
for i = 1:nfms
    hogIntIm(:,:,:,i) = compHoGIntIm(imgs(:,:,:,i));
    hofIntIm(:,:,:,i) = compHoFIntIm(flow.hu(:,:,i), flow.hv(:,:,i));
    [mbhxIntIm(:,:,:,i) mbhyIntIm(:,:,:,i)] = compMBHIntIm(flow.u(:,:,i), flow.v(:,:,i));
end

% Now compute the features for each actionlet
dim = params.num_xcells * params.num_ycells * 9;
for i = 1:length(A)
    bbox = A(i).bbox;
    n = size(bbox, 1);
    deltaT = ceil(n / params.num_tcells);
    hog = blockOrientHistDesc2(hogIntIm(:,:,:,A(i).start:A(i).end), [1:n]', ...
        bbox, params.num_xcells, params.num_ycells);
    hof = blockOrientHistDesc2(hofIntIm(:,:,:,A(i).start:A(i).end), [1:n]', ...
        bbox, params.num_xcells, params.num_ycells);
    mbhx = blockOrientHistDesc2(mbhxIntIm(:,:,:,A(i).start:A(i).end), [1:n]', ...
        bbox, params.num_xcells, params.num_ycells);
    mbhy = blockOrientHistDesc2(mbhyIntIm(:,:,:,A(i).start:A(i).end), [1:n]', ...
        bbox, params.num_xcells, params.num_ycells);
    hog1 = zeros(1, dim * params.num_tcells);
    hof1 = hog1; mbhx1 = hog1; mbhy1 = hog1;
    for j = 1:params.num_tcells
        grFlg = (j-1)*dim +1 : j * dim;
        frFlg = (j-1)*deltaT + 1 : min(n, j * deltaT);
        hog1(grFlg) = sum(hog(frFlg, :)); hof1(grFlg) = sum(hof(frFlg, :));
        mbhx1(grFlg) = sum(mbhx(frFlg, :)); mbhy1(grFlg) = sum(mbhy(frFlg, :));
    end
    F(i, :) = [normr(hog1) normr(hof1) normr(mbhx1) normr(mbhy1)];
end

end