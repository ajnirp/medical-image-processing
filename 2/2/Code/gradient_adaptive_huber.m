function output = gradient_adaptive_huber(x_loop, imageNoisy, alpha, gamma)

y1_in = x_loop-circshift(x_loop,1,1);
y1_out = zeros(size(y1_in));
y1_out(abs(y1_in)<=gamma) = y1_in(abs(y1_in)<=gamma);
y1_out(abs(y1_in)>gamma) = gamma*sign(y1_in(abs(y1_in)>gamma));

y2_in = x_loop-circshift(x_loop,-1,1);
y2_out = zeros(size(y2_in));
y2_out(abs(y2_in)<=gamma) = y2_in(abs(y2_in)<=gamma);
y2_out(abs(y2_in)>gamma) = gamma*sign(y2_in(abs(y2_in)>gamma));

y3_in = x_loop-circshift(x_loop,1,2);
y3_out = zeros(size(y3_in));
y3_out(abs(y3_in)<=gamma) = y3_in(abs(y3_in)<=gamma);
y3_out(abs(y3_in)>gamma) = gamma*sign(y3_in(abs(y3_in)>gamma));

y4_in = x_loop-circshift(x_loop,-1,2);
y4_out = zeros(size(y4_in));
y4_out(abs(y4_in)<=gamma) = y4_in(abs(y4_in)<=gamma);
y4_out(abs(y4_in)>gamma) = gamma*sign(y4_in(abs(y4_in)>gamma));

beta = y1_out + y2_out + y3_out + y4_out;
output = (1-alpha)*2*(x_loop-imageNoisy) + alpha*beta;