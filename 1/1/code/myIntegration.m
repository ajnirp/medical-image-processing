% ------Questions to still answer--------
%%
%
% * How do we find out the delta s
% * Find out a better method to get limits of s. The current method is
% highly inefficient(but maybe we could still manage since its highly
% parallel
% * What would be the ideal image interpolation scheme?(I think linear
% should be fine since our image is a phantom and the input data itself has
% only piecewise constant parts in it
% 

%% Line Integration
% Function performs integration of the image intensities along lines Lt,θ
% that are parametrized by t and θ, where t ∈ (−∞, ∞) is the (signed)
% distance of the line from the origin and θ ∈ [0, 180) is the angle (in
% degrees) made by the upward-pointing normal to the line with the +X axis.
function integration_output = myIntegration(inputImage, t, theta, s_delta)
[M, N] = size(inputImage);

% We know that in Matlab, the minimum value of x and y that the image
% coordinates can take is 1. Also the maximum value that the x coordinate
% can take is M, the maximum value that the y coordinate can take is N.

s1 = (t*cos(theta) - 1)/sin(theta);
s2 = (t*cos(theta) - M)/sin(theta);
s3 = (1-t*sin(theta))/cos(theta);
s4 = (N-t*sin(theta))/cos(theta);
%%
% $x(s) = t*cos(theta)-s*sin(theta);$
%
% $y(s) = t*sin(theta)+s*cos(theta);$

% ----Get these values somehow------
% s_delta = 0.1;

sinitial = min([s1, s2, s3, s4]);
sinitial = max([sinitial,-max(M,N)/s_delta]);
sfinal = max([s1, s2, s3, s4]);
sfinal = min([sfinal,max(M,N)/s_delta]);
s = sinitial:s_delta:sfinal;
n = 1:length(s);

x(n) = M/2 + t*cos(theta)-s*sin(theta);
y(n) = N/2 + t*sin(theta)+s*cos(theta);

[X,Y] = meshgrid(1:M, 1:N);
interpolated_values = interp2(X, Y, inputImage, x(n), y(n), 'linear', 0);
integration_output = sum(interpolated_values) * s_delta;