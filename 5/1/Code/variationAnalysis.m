function [V,D] = variationAnalysis(transposed)
[numOfPoints, ~, numOfPointSets] = size(transposed);
flattened = reshape(transposed, [2*numOfPoints numOfPointSets]);
% temp = shifted(:,:,174); temp = temp'; temp = temp(:)';
% assert(isequal(flattened(:,174)', temp));
shapeCov = zeros([2*numOfPoints 2*numOfPoints]);
for i = 1 : numOfPointSets
   temp = flattened(:,i);
   shapeCov = shapeCov + (temp * temp');
end
shapeCov = shapeCov / numOfPointSets;
[V,D] = eig(shapeCov);