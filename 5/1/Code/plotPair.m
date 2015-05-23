function plotPair(p1,p2)
temp1 = repmat([1;10], [1 32])';
temp2 = [p1,p2];
scatter(temp2(1,:), temp2(2,:), [], temp1(:));