%% Question 2b
S0 = phantom(256);
mask = fspecial ('gaussian', 11, 1);
S1 = conv2 (S0, mask, 'same');
mask = fspecial ('gaussian', 51, 5);
S5 = conv2 (S0, mask, 'same');
% Note S1 is blurred with Gaussian of sigma 1, S5 is blurred with sigma 5

figure();
subplot(1,3,1), imshow(S0), title('Input image')
subplot(1,3,2), imshow(S1), title('Smoothed with gaussian sigma 1')
subplot(1,3,3), imshow(S5), title('Smoothed with gaussian sigma 5')

input_array = [S0, S1, S5];
theta = 0:3:177;
for i = 1:3
    if (i==1)
        image = S0;
        title_name = 'Input Image Filter operation';
        a_number = 0;
    elseif (i==2)
        image = S1;
        title_name = 'S1 smoothed Image Filter operation';
        a_number = 1;
    elseif (i==3)
        image = S5;
        title_name = 'S5 smoothed Image Filter operation';
        a_number = 5;
    end
    
    [A] = radon(image, theta);
    myFilter(1,  
%     extracted_image = iradon(A, theta, 'linear', 'Ram-Lak', 1.0, size(image,1));
    figure(); imshow(extracted_image); title(title_name)
    % RRMSE calculation
    RRMSE = sqrt(sum(sum((extracted_image - image).^2)))/sqrt(sum(sum(image.^2)));
    fprintf('The value of RRMSE for image S%d is %d \n', a_number, RRMSE);
end

%% Observation
% 
% * We observe that RRMSE for image S0>S1>S5
% * Basically, as the image gets more smoothed and smoothed, the RRMSE is
% decreased.
% * This is exactly as we would predict since as we smooth the image with a
% Gaussian filter, the edges in the image become weak and the noise content
% also goes down.
% * Because of these 2 reasons, the inverse radon performs much better.
% This is because the Ram Lak Filter is a linear filter and cuts off
% frequencies higher than a certain threshold. So if the input image does
% not have those higher frequency components, then the reconstructed error
% will be less. Noise and sharp edges introduce high frequencies. Hence the
% error was maximum for image S0.  But since image S5 is the smoothest, the
% error for it is the least.
% 
