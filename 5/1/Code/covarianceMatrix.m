function output = covarianceMatrix(psets)
meanMatrix = estimateMean(psets);
meanReshaped = reshape(meanMatrix, [1 1 300]);
shifted = psets - meanReshaped;
output = shifted