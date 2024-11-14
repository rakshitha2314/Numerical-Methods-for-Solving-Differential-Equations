function [x, y] = linear_finite_difference(p, q, r, a, b, alpha, beta, N)
    % LINEAR_FINITE_DIFFERENCE Solves a linear BVP using the finite difference method
    % p, q, r: Function handles for p(x), q(x), and r(x) in the differential equation
    %          y'' = p(x) y' + q(x) y + r(x)
    % a, b:    The interval [a, b] where the BVP is defined
    % alpha:   Boundary condition y(a) = alpha
    % beta:    Boundary condition y(b) = beta
    % N:       Number of subintervals (number of interior points is N)
    
    h = (b - a) / (N + 1);  % Step size
    x = linspace(a, b, N+2); % Grid points, including boundary points
    
    % Initialize the system matrix A and vector B
    A = zeros(N, N);
    B = zeros(N, 1);
    
    % Loop over each interior point and set up the finite difference equations
    for i = 1:N
        xi = a + i * h; % Current x value at interior point
        
        % Fill in the matrix A for the finite difference scheme
        A(i, i) = -2 + h^2 * q(xi); % Diagonal element
        if i > 1
            A(i, i-1) = 1 - (h/2) * p(xi); % Lower diagonal element
        end
        if i < N
            A(i, i+1) = 1 + (h/2) * p(xi); % Upper diagonal element
        end
        
        % Fill in the vector B with boundary conditions
        B(i) = h^2 * r(xi);
        if i == 1
            B(i) = B(i) - (1 - (h/2) * p(xi)) * alpha; % Adjust for y(a) = alpha
        elseif i == N
            B(i) = B(i) - (1 + (h/2) * p(xi)) * beta; % Adjust for y(b) = beta
        end
    end
    
    % Solve the linear system A * y_interior = B
    y_interior = A \ B;
    
    % Construct the full solution, including boundary values
    y = [alpha; y_interior; beta];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters and exact solutions for both questions

% Question 8: Define BVP and exact solution
p8 = @(x) 0;
q8 = @(x) -4;
r8 = @(x) 4 * x;
a8 = 0;
b8 = 1;
alpha8 = 0;
beta8 = 2;
y_exact8 = @(x) x + (exp(2) * (exp(2 * x) - exp(-2 * x))) / (exp(4) - 1);

% Question 9: Define BVP and exact solution
p9 = @(x) 0;
q9 = @(x) -100;
r9 = @(x) 0;
a9 = 0;
b9 = 1;
alpha9 = 1;
beta9 = exp(-10);
y_exact9 = @(x) exp(-10 * x);

% Helper function to plot results for each question
function plot_combined_results(x_a, y_a, h_a, x_b, y_b, h_b, y_exact, question_num)
    x_exact = linspace(min(x_a), max(x_a), 100);
    y_exact_values = y_exact(x_exact);
    
    figure;
    plot(x_a, y_a, 'bo-', 'DisplayName', sprintf('Approximate Solution, h = %.2f', h_a), LineWidth=1);
    hold on;
    plot(x_b, y_b, 'rx-', 'DisplayName', sprintf('Approximate Solution, h = %.2f', h_b), LineWidth=1);
    plot(x_exact, y_exact_values, 'k--', 'LineWidth', 1, 'DisplayName', 'Exact Solution');
    
    xlabel('x');
    ylabel('y');
    legend;
    title(sprintf('Finite Difference Solution vs Exact Solution (Question %d)', question_num));
    hold off;
end

% Question 8: Apply Linear Finite Difference Method with h = 1/2 and h = 1/4
fprintf('Question 8: Solving with h = 1/2\n');
[x8a, y8a] = linear_finite_difference(p8, q8, r8, a8, b8, alpha8, beta8, 1); % h = 1/2

fprintf('Question 8: Solving with h = 1/4\n');
[x8b, y8b] = linear_finite_difference(p8, q8, r8, a8, b8, alpha8, beta8, 3); % h = 1/4

% Plot both h = 1/2 and h = 1/4 on the same graph for Question 8
plot_combined_results(x8a, y8a, 1/2, x8b, y8b, 1/4, y_exact8, 8);

% Question 9: Apply Linear Finite Difference Method with h = 1/10 and h = 1/20
fprintf('Question 9: Solving with h = 1/10\n');
[x9a, y9a] = linear_finite_difference(p9, q9, r9, a9, b9, alpha9, beta9, 9); % h = 1/10

fprintf('Question 9: Solving with h = 1/20\n');
[x9b, y9b] = linear_finite_difference(p9, q9, r9, a9, b9, alpha9, beta9, 19); % h = 1/20

% Plot both h = 1/10 and h = 1/20 on the same graph for Question 9
plot_combined_results(x9a, y9a, 1/10, x9b, y9b, 1/20, y_exact9, 9);