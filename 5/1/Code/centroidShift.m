% Find the centroid of a pointset and subtract it from all points in the
% pointset
function output = centroidShift(ps)
numPoints = size(ps, 2);
centroid = sum(ps, 2) / numPoints;
output = ps - repmat(centroid, [1 numPoints]);