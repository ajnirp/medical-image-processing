% Plot all the pointsets with the assigned colors
function plotAllSets(pointSets, filename, plotTitle)

colours = makeColours(pointSets);

numOfPoints = size(pointSets, 2);
numOfPointSets = size(pointSets, 3);

reshaped = reshape(pointSets, [2 numOfPoints*numOfPointSets]);
% assert(isequal(reshaped(:,32*2+1:32*3), pointSets(:,:,3)));

title(plotTitle);
h = scatter(reshaped(1,:,:), reshaped(2,:,:), [], colours, 'o');
saveas(h, strcat('../Images/', filename), 'png');