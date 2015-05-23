%% Denoising a Phantom Magnetic Resonance Image
% We want to denoise the image using Bayesian image-denoising algorithm
% that uses a noise model coupled with a MRF prior that uses a 4-neighbor
% neighborhood system
clear all
load('assignmentImageDenoisingPhantom.mat');
alpha = 0.5;
tow = 0.1;
gamma = 1;
x_old = imageNoisy;
iteration_count = 1;
alpha = 0.05:0.05:0.95;
find_optimal_alpha_matrix = [inf, inf, inf];
optim_fn_iteration_variation = [inf, inf];
for gamma = 0.01
    for alpha = 0.9
        tow = 0.1;
        x_old = imageNoisy;
        iteration_count = 1;
        while(tow>0.01)
            g_x_old = gradient_adaptive_huber(x_old, imageNoisy, alpha, gamma);
            x = x_old - tow*g_x_old;
            g_x_new = gradient_adaptive_huber(x, imageNoisy, alpha, gamma);
            if(sum(sum(abs(g_x_old))) > (sum(sum(abs(g_x_new)))))
                tow = 1.1*tow;
                x_old = x;
            else
                tow = 0.5*tow;
                x = x_old;
            end
            optim_fn_iteration_variation = [optim_fn_iteration_variation; [iteration_count, sum(sum(abs(g_x_new)))]];
            iteration_count = iteration_count + 1;
        end
        find_optimal_alpha_matrix = [find_optimal_alpha_matrix; [alpha,gamma, RRMSE(imageNoiseless, x)]];
        fprintf('Number of iterations the gradient descent ran was %f \n', iteration_count);
        RRMSE(imageNoiseless, x)
    end
end
imshow(abs(x))
%%
% We can verify from above that the optimal value of alpha is 0.9, gamma is
% 0.01, for RRMSE equal to 0.0662
%
% For alpha = 0.9, gamma = 0.01, RRMSE = 0.0662
%
% For alpha = 1, gamma = 0.01, RRMSE = 0.8112
%
% For alpha = 0.72, gamma = 0.01, RRMSE = 0.1949
%
% For alpha = 0.9, gamma = 0.012, RRMSE = 0.0703
%
% For alpha = 0.9, gamma = 0.008, RRMSE = 0.0787
%