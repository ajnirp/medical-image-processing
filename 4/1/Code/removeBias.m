function output = removeBias(memb, class)
[r c] = size(memb(:,:,1));
output = zeros([r c]);
for i = 1:r
    for j = 1:c
        output(i,j) = dot(reshape(memb(i,j,:), [3 1]), class);
    end
end