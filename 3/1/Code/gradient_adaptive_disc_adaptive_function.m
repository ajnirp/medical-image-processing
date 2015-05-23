function output = gradient_adaptive_disc_adaptive_function(x_loop, imageNoisy, alpha, gamma, mask)

y1_in = x_loop-circshift(x_loop,1,1);
y1_out = (gamma*y1_in)./(gamma+abs(y1_in));

y2_in = x_loop-circshift(x_loop,-1,1);
y2_out = (gamma*y2_in)./(gamma+abs(y2_in));

y3_in = x_loop-circshift(x_loop,1,2);
y3_out = (gamma*y3_in)./(gamma+abs(y3_in));


y4_in = x_loop-circshift(x_loop,-1,2);
y4_out = (gamma*y4_in)./(gamma+abs(y4_in));

beta = y1_out + y2_out + y3_out + y4_out;
output = (1-alpha)*2*(ifft2(mask .* fft2(x_loop)) - ifft2(mask .* imageNoisy)) + alpha*beta;