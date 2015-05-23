% Observe the variation of the shape along the k-th mode
function output = reconstructHelper(meanMatrix, V, D, k, mult)

numOfPoints = size(meanMatrix, 2);
N = 2*numOfPoints;

meanReshaped = reshape(meanMatrix', [1 2*numOfPoints]);
totalShift = zeros(size(meanReshaped));

% for i = 1 : k
for i = k : k
    j = N - i + 1;
    shift = V(:,j)';
    coeff = mult*sqrt(D(j,j));
    totalShift = totalShift + coeff*shift;
end

output = meanReshaped + totalShift;
output = [output(1:numOfPoints) ; output(numOfPoints+1:2*numOfPoints)];