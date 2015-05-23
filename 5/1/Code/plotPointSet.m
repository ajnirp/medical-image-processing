function plotPointSet(ps, ti, fn)
h = scatter(ps(1,:), ps(2,:));
title(ti);
saveas(h, strcat('../Images/', fn), 'png');