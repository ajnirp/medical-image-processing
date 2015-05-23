function output_image = ApplyFilter(input_image, filter)

theta = 0:3:177;
A = radon(input_image, theta);

ft_R = fftshift(fft(A,[],1),1);
filteredProj = ft_R .* filter;
filteredProj = ifftshift(filteredProj,1);
ift_R = real(ifft(filteredProj,[],1));
output_image = iradon(ift_R, theta, 'linear', 'none', 1.0, size(input_image,1));