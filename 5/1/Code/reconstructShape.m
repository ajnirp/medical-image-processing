% Observe the variation of the shape along the k-th mode
function output = reconstructShape(meanMatrix, V, D, k)
numOfPoints = size(meanMatrix, 2);
output = zeros([2 numOfPoints 5]);
mults = [0 1 -1 2 -2];
for i = 1 : 5
    output(:,:,i) = reconstructHelper(meanMatrix, V, D, k, mults(i));
end