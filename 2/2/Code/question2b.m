%% Denoising a Phantom Magnetic Resonance Image
% We want to denoise the image using Bayesian image-denoising algorithm
% that uses a noise model coupled with a MRF prior that uses a 4-neighbor
% neighborhood system
clear all
load('assignmentImageDenoisingBrainNoisy.mat');
alpha = 0.5;
tow = 0.1;
x_old = imageNoisy;
iteration_count = 1;
alpha = 0.05:0.05:0.95;
find_optimal_alpha_matrix = [inf, inf];
optim_fn_iteration_variation = [inf, inf];
for alpha = 0.3
    tow = 0.1;
    x_old = imageNoisy;
    iteration_count = 1;
    while(tow>0.01)
        g_x_old = gradient_quadfunction(x_old, imageNoisy, alpha);
        x = x_old - tow*g_x_old;
        g_x_new = gradient_quadfunction(x, imageNoisy, alpha);
        if(sum(sum((abs(g_x_old)>abs(g_x_new)))) > sum(sum((abs(g_x_old)<abs(g_x_new)))))
            tow = 1.1*tow;
            x_old = x;
        else
            tow = 0.5*tow;
            x = x_old;
        end
        optim_fn_iteration_variation = [optim_fn_iteration_variation; [iteration_count, sum(sum(abs(g_x_new)))]];
        iteration_count = iteration_count + 1;
    end
    fprintf('Number of iterations the gradient descent ran was %f \n', iteration_count);
    
    figure();
    imshow(abs(x))
    s = 'Optimal Quadratic Prior Image for ALPHA =  ';
    s = strcat(s,num2str(alpha));
    title(s);
end

%%
% We can verify from visual inspection that the optimal value of alpha is 0.3

imshow(abs(x));