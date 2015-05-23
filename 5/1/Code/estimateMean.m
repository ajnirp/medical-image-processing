function output = estimateMean(psets)
avg = sum(psets, 3) / size(psets, 3);
% avg = sum(psets, 3) / numOfPointSets;
output = avg / norm(avg, 2);