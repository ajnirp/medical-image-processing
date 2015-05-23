function output = finalClassEstimate(memb, class)
[r c] = size(memb(:,:,1));
output = zeros([r c]);
for i = 1:r
    for j = 1:c
        membs = memb(i,j,:);
        assigned = find(membs == max(membs));
        output(i,j) = class(assigned);
    end
end