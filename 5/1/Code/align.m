% Code to align two pointsets of equal cardinality via similarity
% transformations
function res = align(ps1, ps2)
% shifting
sh1 = centroidShift(ps1);
sh2 = centroidShift(ps2);
% scaling
sc1 = sh1 / norm(sh1, 2);
sc2 = sh2 / norm(sh2, 2);
% we assume all points are equally important i.e. weight matrix W = I
[U,~,V] = svd(sc1 * sc2');
rot = V*U';
if det(rot) - 1 > 1e-5
    display('need to negate');
    rot = V*[1 0; 0 -1]*U';
end
res = rot*sc1;