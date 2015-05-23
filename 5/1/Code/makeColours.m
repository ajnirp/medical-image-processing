function colours = makeColours(pointSets)

numOfPointSets = size(pointSets, 3);
numOfPoints = size(pointSets, 2);

colours = 1 : 32 : 32*numOfPointSets;
colours = repmat(colours, [numOfPoints 1]);

colours = colours(:)';