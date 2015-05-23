%% Denoising a Phantom Magnetic Resonance Image
% We want to denoise the image using Bayesian image-denoising algorithm
% that uses a noise model coupled with a MRF prior that uses a 4-neighbor
% neighborhood system
clear all
load('assignmentImageDenoisingBrainNoisy.mat');
alpha = 0.5;
tow = 0.1;
gamma = 1;
x_old = imageNoisy;
iteration_count = 1;
alpha = 0.05:0.05:0.95;
find_optimal_alpha_matrix = [inf, inf, inf];
optim_fn_iteration_variation = [inf, inf];
for alpha = 0.6:0.1:0.9
    for gamma = 0.015
        tow = 5;
        x_old = imageNoisy;
        iteration_count = 1;
        RRMSE_new = 10;
        RRMSE_old = 5;
        while(tow>0.01 && abs(RRMSE_new - RRMSE_old)>(0.0001*RRMSE_new))
            g_x_old = gradient_adaptive_disc_adaptive_function(x_old, imageNoisy, alpha, gamma);
            x = x_old - tow*g_x_old;
            g_x_new = gradient_adaptive_disc_adaptive_function(x, imageNoisy, alpha, gamma);
            if(sum(sum(abs(g_x_old))) > (sum(sum(abs(g_x_new)))))
                tow = 1.1*tow;
                x_old = x;
            else
                tow = 0.5*tow;
                x_old = x;
            end
            optim_fn_iteration_variation = [optim_fn_iteration_variation; [iteration_count, sum(sum(abs(g_x_new)))]];
            %             figure(1);
            iteration_count = iteration_count + 1;
%             RRMSE_new = RRMSE(imageNoiseless, x);
        end
        %         find_optimal_alpha_matrix = [find_optimal_alpha_matrix; [alpha,gamma, RRMSE(imageNoiseless, x)]];
        %         fprintf('Number of iterations the gradient descent ran was %f \n', iteration_count);
        %                 RRMSE(imageNoiseless, x)
    end
    s = 'Optimal Adaptive function Prior Image for ALPHA =  ';
    s = strcat(s,num2str(alpha),'and gamma = ');
    s = strcat(s, num2str(gamma));
    figure();
    imshow(abs(x))
    title(s);
    pause();
end
%%
% We can verify from above that the optimal value of alpha is 0.9, gamma is
% 0.015, for RRMSE equal to 0.0854
%
% For alpha = 0.9, gamma = 0.015, RRMSE = 0.0854
%
% For alpha = 0.99, gamma = 0.015, RRMSE = 0.3767
%
% For alpha = 0.72, gamma = 0.015, RRMSE = 0.1775
%
% For alpha = 0.9, gamma = 0.018, RRMSE = 0.0873
%
% For alpha = 0.9, gamma = 0.012, RRMSE = 0.0875
%