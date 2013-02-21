% explode the frame and blur it
% @input: input frame that is single channel
% @quant_num: number of quantization
% @sigma_blur: std for gaussian blur

function f = getFrameExpBlur(input,quant_num,sigma_blur)

if size(input,3)>1
    error('input frame must have a single channel!');
end

% explode
step = 256/quant_num;
f=zeros(size(input,1),size(input,2),quant_num);
thresh=0;
for i=1:quant_num
    f(:,:,i)=(input>=thresh & input<thresh+step);
    thresh=thresh+step;
end

% blur

h=fspecial('gaussian',[min(20,sigma_blur),min(20,sigma_blur)],sigma_blur);
for i=1:quant_num
    f(:,:,i) = imfilter(f(:,:,i),h,'same');
end

