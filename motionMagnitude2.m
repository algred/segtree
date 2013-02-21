function M = motionMagnitude2(flow)

M = zeros(size(flow.hu));
for i = 1:size(flow.hu, 3)
    M(:,:,i) = sqrt(flow.u(:,:,i).*flow.u(:,:,i) + flow.v(:,:,i).*flow.v(:,:,i));
end

end