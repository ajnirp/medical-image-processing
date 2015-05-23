%% X Ray computed Tomography : Radon Transform

%%
% We use the Shepp-Logan phantom image here
Fxy = phantom(128);
imshow(Fxy);
title('Shepp-Logan Phantom image')

%%
% We will have to shift the coordinate frame to the center since the
% default matlab convention places the origin at one of the edges of the
% image.
