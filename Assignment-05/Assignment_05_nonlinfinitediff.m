function [x, y] = nonlinear_finite_difference(f, a, b, alpha, beta, N, tol, max_iter)
    % NONLINEAR_FINITE_DIFFERENCE Solves a nonlinear BVP using the finite difference method
    % f: Function handle for the nonlinear term f(x, y, y')
    % a, b: The interval [a, b] where the BVP is defined
    % alpha, beta: Boundary conditions y(a) = alpha, y(b) = beta
    % N: Number of subintervals (number of interior points is N)
    % tol: Tolerance for convergence in Newton's method
    % max_iter: Maximum number of iterations for Newton's method
    
    h = (b - a) / (N + 1);  % Step size
    x = linspace(a, b, N+2); % Grid points, including boundary points
    
    % Initial guess for y (start with linear interpolation between boundaries)
    y = linspace(alpha, beta, N+2)';  % Include boundary conditions in y
    
    % Define function to compute residuals for Newton's method
    function R = residual(y)
        % Set up the residual vector R
        R = zeros(N, 1);
        for i = 2:N+1
            y_prev = y(i-1);
            y_curr = y(i);
            y_next = y(i+1);
            
            % Approximations for y'' and y' at x_i
            y_double_prime = (y_next - 2*y_curr + y_prev) / h^2;
            y_prime = (y_next - y_prev) / (2 * h);
            
            % Nonlinear equation at each interior point
            R(i-1) = y_double_prime - f(x(i), y_curr, y_prime);
        end
    end

    % Newton's method iteration
    for iter = 1:max_iter
        % Compute the residual vector
        R = residual(y);
        
        % Check for convergence
        if norm(R, Inf) < tol
            disp(['Converged in ', num2str(iter), ' iterations']);
            return;
        end
        
        % Approximate the Jacobian matrix J
        J = zeros(N, N);
        delta = 1e-8; % Small perturbation for finite difference
        for j = 1:N
            y_perturbed = y;
            y_perturbed(j+1) = y(j+1) + delta; % Perturb y_j
            R_perturbed = residual(y_perturbed);
            J(:, j) = (R_perturbed - R) / delta;
        end
        
        % Update step: Solve J * delta_y = -R
        delta_y = J \ -R;
        
        % Update interior values of y
        y(2:N+1) = y(2:N+1) + delta_y;
    end
    
    error('Newton method failed to converge within maximum iterations');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Helper function to plot results for each question
function plot_nonlinear_results(x, y_approx, y_exact, question_num)
    x_exact = linspace(min(x), max(x), 100);
    y_exact_values = y_exact(x_exact);
    
    figure;
    plot(x, y_approx, 'bo-', 'DisplayName', 'Approximate Solution', LineWidth=1);
    hold on;
    plot(x_exact, y_exact_values, 'r--', 'LineWidth', 1, 'DisplayName', 'Exact Solution');
    
    xlabel('x');
    ylabel('y');
    legend;
    title(sprintf('Nonlinear Finite Difference Solution vs Exact Solution (Question %d)', question_num));
    hold off;
end

% Question 10
fprintf('Question 10: Solving with h = 0.5\n');
f10 = @(x, y, yp) -(yp^2) - y + log(x);
a10 = 1;
b10 = 2;
alpha10 = 0;
beta10 = log(2);
N10 = (b10 - a10) / 0.5 - 1; % Calculate N from h = 0.5
tol10 = 1e-6; % Tolerance for Newton's method
max_iter10 = 20; % Maximum number of iterations
y_exact10 = @(x) log(x);

[x10, y10] = nonlinear_finite_difference(f10, a10, b10, alpha10, beta10, N10, tol10, max_iter10);
plot_nonlinear_results(x10, y10, y_exact10, 10);

% Question 11
fprintf('Question 11: Solving with N = 10, TOL = 10^(-4)\n');
f11 = @(x, y, yp) -exp(-2 * y);
a11 = 1;
b11 = 2;
alpha11 = 0;
beta11 = log(2);
N11 = 9; % N is specified as 10
tol11 = 1e-4; % Tolerance for Newton's method
max_iter11 = 20; % Maximum number of iterations
y_exact11 = @(x) log(x);

[x11, y11] = nonlinear_finite_difference(f11, a11, b11, alpha11, beta11, N11, tol11, max_iter11);
plot_nonlinear_results(x11, y11, y_exact11, 11);

% Question 12
fprintf('Question 12: Solving with N = 10\n');
f12 = @(x, y, yp) yp * cos(x) - y * log(y);
a12 = 0;
b12 = pi / 2;
alpha12 = 1;
beta12 = exp(1);
N12 = 9; % N is specified as 10
tol12 = 1e-6; % Tolerance for Newton's method
max_iter12 = 20; % Maximum number of iterations
y_exact12 = @(x) exp(sin(x));

[x12, y12] = nonlinear_finite_difference(f12, a12, b12, alpha12, beta12, N12, tol12, max_iter12);
plot_nonlinear_results(x12, y12, y_exact12, 12);