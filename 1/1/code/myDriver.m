function integration_result = myRadonTrans()
% In the output integration_result is a matrix of size 37*36, since t takes
% 37 values and theta takes 36. the value of integration_result in i'th row
% and j'th column corresponds to the value of radon transform with t = i'th
% value from the sequence of -90:5:90 and theta = j'th value from the
% sequence 0:5:175(these values were as asked in the question)
%%
% We use the Shepp-Logan phantom image here
Fxy = phantom(128);
imshow(Fxy);
title('Shepp-Logan Phantom image')

t = -90:5:90;
theta = 0:5:175;
integration = zeros(length(t),length(theta));

for i = 1:length(t)
    for j = 1:length(theta)
        integration(i,j) = myIntegration(Fxy, t(i), theta(j), 1);
    end
end
integration_result = integration;