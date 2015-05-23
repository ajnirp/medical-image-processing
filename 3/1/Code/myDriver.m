%% Reconstructing a Phantom Magnetic Resonance Image - Quadratic Prior
% We want to reconstruct the image using Bayesian image-denoising algorithm
% that uses a noise model coupled with a MRF prior that uses a 4-neighbor
% neighborhood system
clear all
load('assignmentImageReconstructionPhantom.mat');
RRMSE(ifft2(imageKspaceData), imageNoiseless)
find_optimal_alpha_matrix = [inf, inf];
optim_fn_iteration_variation = [inf, inf];
for alpha = 0.33;
    x_old = ifft2(imageKspaceData);
    tow = 0.1;
    iteration_count = 1;
    while(tow>0.005)
        g_x_old = gradient_quadfunction(x_old, imageKspaceData, alpha, imageKspaceMask);
        fourier_term = 2*(x_old - ifft2(imageKspaceData));
        x = x_old - tow*fourier_term - tow*g_x_old;
        g_x_new = gradient_quadfunction(x, imageKspaceData, alpha, imageKspaceMask);
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
    imwrite(abs(x), strcat('../Images/quad-prior-', num2str(alpha), '.png'), 'png');
    find_optimal_alpha_matrix = [find_optimal_alpha_matrix; [alpha, RRMSE(imageNoiseless, x)]];
    fprintf('Number of iterations the gradient descent ran was %f \n', iteration_count);
end
% imwrite(ifft2(imageKspaceData), '../Images/original-data.tiff', 'tiff');
quad_reconstructed = x;
% figure();
% imshow(abs(x))
% title('Best reconstructed image using quadratic prior')
%%
% We can verify from above that the optimal value of alpha is 0.33, for
% RRMSE equal to 0.2435
%
% For alpha = 0.33, RRMSE = 0.2435
% For alpha = 0.396, RRMSE = 0.2441
% For alpha = 0.264, RRMSE = 0.2440
%%
h=figure();
plot(optim_fn_iteration_variation(2:end,1), optim_fn_iteration_variation(2:end,2));
title('Variation of optimality criteria quadratic prior(derivative of MAP)')
xlabel('Number of iterations')
ylabel('Value of optimality critera(that we want to minimise')
saveas(h,'../Images/quad-prior-plot','png');
%% Reconstructing a Phantom Magnetic Resonance Image - Huber function
% We want to reconstruct the image using Bayesian image-denoising algorithm
% that uses a noise model coupled with a MRF prior that uses a 4-neighbor
% neighborhood system
% clear all
load('assignmentImageReconstructionPhantom.mat');
iteration_count = 1;
find_optimal_alpha_matrix = [inf, inf, inf];
optim_fn_iteration_variation = [inf, inf];
for gamma = 0.03
    for alpha = 0.73
        tow = 0.1;
        x_old = ifft2(imageKspaceData);
        iteration_count = 1;
        while(tow>0.01)
            g_x_old = gradient_adaptive_huber(x_old, imageKspaceData, alpha, gamma, imageKspaceMask);
            x = x_old - tow*g_x_old;
            g_x_new = gradient_adaptive_huber(x, imageKspaceData, alpha, gamma, imageKspaceMask);
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
        imwrite(abs(x), strcat('../Images/huber-prior-', num2str(gamma), num2str(alpha), '.png'), 'png');
        find_optimal_alpha_matrix = [find_optimal_alpha_matrix; [alpha,gamma, RRMSE(imageNoiseless, x)]];
        fprintf('Number of iterations the gradient descent ran was %f \n', iteration_count);
    end
end
% figure();
huber_reconstructed = x;
% imshow(abs(x))
% title('Best reconstructed image using huber prior')
%%
h=figure();
plot(optim_fn_iteration_variation(2:end,1), optim_fn_iteration_variation(2:end,2));
title('Variation of optimality criteria huber prior(derivative of MAP)')
xlabel('Number of iterations')
ylabel('Value of optimality critera(that we want to minimise')
saveas(h,'../Images/huber-prior-plot','png');
%%
% We can verify from above that the optimal value of alpha is 0.73, gamma is
% 0.03, for RRMSE equal to 0.1988
%
% For alpha = 0.73, gamma = 0.03, RRMSE = 0.1965
%
% For alpha = 0.876, gamma = 0.03, RRMSE = 0.1983
% For alpha = 0.584, gamma = 0.03, RRMSE = 0.2081
%
% For alpha = 0.73, gamma = 0.036, RRMSE = 0.1966
% For alpha = 0.73, gamma = 0.024, RRMSE = 0.1974
%
%% Reconstructing a Phantom Magnetic Resonance Image - Discontinuity Adaptive Function
% We want to reconstruct the image using Bayesian image-denoising algorithm
% that uses a noise model coupled with a MRF prior that uses a 4-neighbor
% neighborhood system
% clear all
load('assignmentImageReconstructionPhantom.mat');
find_optimal_alpha_matrix = [inf, inf, inf];
optim_fn_iteration_variation = [inf, inf];
for gamma = 0.00144
    for alpha = 0.999
        tow = 5;
        x_old = ifft2(imageKspaceData);
        iteration_count = 1;
        RRMSE_new = 10;
        RRMSE_old = 5;
        while(tow>0.01 && abs(RRMSE_new - RRMSE_old)>(0.0001*RRMSE_new))
            g_x_old = gradient_adaptive_disc_adaptive_function(x_old, imageKspaceData, alpha, gamma, imageKspaceMask);
            x = x_old - tow*g_x_old;
            g_x_new = gradient_adaptive_disc_adaptive_function(x, imageKspaceData, alpha, gamma, imageKspaceMask);
            RRMSE_old = RRMSE(imageNoiseless, x_old);
            if(sum(sum(abs(g_x_old))) > (sum(sum(abs(g_x_new)))))
                tow = 1.1*tow;
                x_old = x;
            else
                tow = 0.5*tow;
                x_old = x;
            end
            optim_fn_iteration_variation = [optim_fn_iteration_variation; [iteration_count, sum(sum(abs(g_x_new)))]];
            iteration_count = iteration_count + 1;
            RRMSE_new = RRMSE(imageNoiseless, x);
        end
        imwrite(abs(x), strcat('../Images/adap-prior-', num2str(gamma), num2str(alpha), '.png'), 'png');
        find_optimal_alpha_matrix = [find_optimal_alpha_matrix; [alpha,gamma, RRMSE(imageNoiseless, x)]];
        fprintf('Number of iterations the gradient descent ran was %f \n', iteration_count);
    end
end
adap_reconstructed = x;
% figure();
% imshow(abs(x))
% title('Best reconstructed image using discontinuity adaptive prior')
%%
h=figure();
plot(optim_fn_iteration_variation(2:end,1), optim_fn_iteration_variation(2:end,2));
title('Variation of optimality criteria discontinuity adaptive prior(derivative of MAP)')
xlabel('Number of iterations')
ylabel('Value of optimality critera(that we want to minimise')
saveas(h,'../Images/adap-prior-plot','png');
%%
% We can verify from above that the optimal value of alpha is 0.995, gamma is
% 0.0012, for RRMSE equal to 0.1863
%
% For alpha = 0.995, gamma = 0.0012, RRMSE = 0.1863
%
% For alpha = 0.796, gamma = 0.0012, RRMSE = 0.2513
% For alpha = 0.999, gamma = 0.0012, RRMSE = 0.1865
%
% For alpha = 0.995, gamma = 0.00096, RRMSE = 0.1863
% For alpha = 0.995, gamma = 0.00144, RRMSE = 0.1866
%%
h=figure('units','normalized','outerposition',[0 0 1 1]);
title('All images compared')
subplot(1,5,1), imshow(ifft2(imageKspaceData)), title('Noisy Image')
subplot(1,5,2), imshow(abs(quad_reconstructed)), title('Reconstructed using Quadratic prior')
subplot(1,5,3), imshow(abs(huber_reconstructed)), title('Reconstructed using Huber prior')
subplot(1,5,4), imshow(abs(adap_reconstructed)), title('Reconstructed using Adaptive Function prior')
subplot(1,5,5), imshow(abs(imageNoiseless)), title('Noiseless Image'), colorbar
saveas(h,'../Images/combined','png');