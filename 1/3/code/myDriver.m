% input values
g = [1 0 ; 0.866 0.5 ; 0.5 0.866 ; 0 1 ; -0.5 0.866 ; -0.866 0.5];
g = g';

% observations
% note: we discard noise by taking the observations to be the magnitude of
% the complex observations
S_actual = [0.5045-0.0217i 0.6874+0.0171i 0.3632+0.1789i 0.3483+0.1385i 0.2606-0.0675i 0.2407+0.1517i];
S = abs(S_actual)';
% the b parameter in the diffusion tensor model
b = 0.1;
% the step size
step = 0.5;

% Let D = LL' where L = [a 0 ; c b]
% i.e. D = [a^2 ac ; ac c^2+b^2]
% Then the objective function is:
% the sum of (i = 1 to 6) of the square
% S(i) - k*exp((a*x + c*y)^2 + (b*y)^2)
% where x = g(i,1), y = g(i,2)

% This is a least squares optimization problem
% And the residual function is non-linear since it involves exp
% So we use the Levenberg-Marqauardt algorithm

% As we can see from the above residual function
% we have three unknowns, a b and c (x and y are elements of g and thus
% constant). So our Jacobian is 3x3
% Let L = [a 0;b c]

% Update step: beta_new = L - s(J'J + lambda*diag(J'J))^-1(J'r)
% where s is a step size (< 1)

% beta = [.1 .1 .1]';

L = rand(2,2); % the lower triangular matrix obtained by the cholesky decomp of D
L(1,2) = 0; % enforcing lower triangularity

iteration = 500; % we run levenberg-marqaudt for a fixed number of iterations
er = 1e5; % current error in the iteration

% arrays in which we will store the quantities to be plotted at the end of
% the script
errors = [];
diag1_values = [];
diag2_values = [];
offdiag_values = [];

for j = 1:iteration
    D = L * L';
    x = diag(g' * D * g);
    t = exp(-b * x); % estimate
    er = S - t; % error = observed - estimate
    derivative = -2 * b * (er .* t);
    Prod = L'*g;
    temp1 = [Prod(1,:) ; Prod(1,:) ; Prod(2,:)];
    temp2 = [g(1,:) ; g(2,:) ; g(2,:)];
    J = (temp1 .* temp2)';
    
    % compute the final step size
    JtrspJ = J' * J;
    delta = step*eye(3)/(JtrspJ + (eye(3).*JtrspJ))*J'*er ;
    
    % update the unknowns i.e. the elements of L
    L = [L(1,1) - delta(1) 0 ; L(2,1) - delta(2) L(2,2) - delta(3)];
    
    % the actual error is the squared norm of 'er'
    er = dot(er, er);
    
    % update D
    D = L * L';
    % projection step to ensure that positive diagonal elements of D do not
    % become negative
    if D(1,1) < 0
        D(1,1) = 0.001;
    end
    if D(2,2) < 0
        D(2,2) = 0.001;
    end
    
    errors(j) = er;
    diag1_values(j) = D(1,1);
    diag2_values(j) = D(2,2);
    offdiag_values(j) = D(2,1);
end

% print the answers to the questions
D = L * L'
[evecs, evals]=eig(D);
evalue_ratio = evals(2,2) / evals(1,1)
evalue_ratio_rounded = round(evalue_ratio)
principal_evector = evecs(:,2)

evecs

% return

% draw all the required plots
figure
subplot(2,2,1)
plot(errors);
title('Errors vs # iterations')
subplot(2,2,2)
plot(diag1_values);
title('D(1,1) vs # iterations')
subplot(2,2,3)
plot(diag2_values);
title('D(2,2) vs # iterations')
subplot(2,2,4)
plot(offdiag_values);
title('D(1,2) == D(2,1) vs # iterations')