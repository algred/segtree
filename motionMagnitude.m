function M = motionMagnitude(flow)

M = zeros(size(flow.hu));
for i = 1:size(flow.hu, 3)
    M(:,:,i) = sqrt(flow.hu(:,:,i).*flow.hu(:,:,i) + flow.hv(:,:,i).*flow.hv(:,:,i));
end

end