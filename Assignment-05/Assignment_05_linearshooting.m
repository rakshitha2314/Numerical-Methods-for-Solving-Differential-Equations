% QUESTION-01

% Define parameters for the BVP
a = 1;          % Initial point
b = 2;          % End point
alpha = 1;      % Boundary condition at x = a
beta = 2;       % Boundary condition at x = b
N = 10;         % Number of steps

% Define the function f(x, y, y') in the form required by the problem
p = @(x) -2 / x;
q = @(x) 2 / x^2;
r = @(x) sin(log(x)) / x;

% Run the Linear Shooting Method
[x_vals, y_vals] = linear_shooting_method(p, q, r, a, b, alpha, beta, N);

% Plot the result
figure;
plot(x_vals, y_vals, 'b-o', LineWidth = 1);
xlabel('x');
ylabel('y');
title('Solution of BVP using Linear Shooting Method (Question1)');
grid on;
hold on;

% Define the exact solution for comparison
c2 = -0.03920701320;
c1 = 1.1392070132;
exact_solution = @(x) c1 * x + c2 / x^2 - (3/10) * sin(log(x)) - (1/10) * cos(log(x));

% Evaluate and plot the exact solution at the same x values
y_exact = arrayfun(exact_solution, x_vals);
plot(x_vals, y_exact, 'r--', LineWidth = 1);
legend('Numerical Solution', 'Exact Solution');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%QUESTION-02
% Define parameters for the BVP
a = 0;  % x = 0
b = 1;  % x = 1
alpha = 0;  % y(0) = 0
beta = 2;  % y(1) = 2
N_half = 2; % Number of steps for h = 1/2 (N = (b - a)/h = 1/2)
N_quarter = 4; % Number of steps for h = 1/4 (N = (b - a)/h = 1/4)

% Define p, q, and r for the BVP y'' = 4(y - x)
p = @(x) 0;  % p(x) = 0
q = @(x) 4;  % q(x) = 4
r = @(x) -4*x;  % r(x) = -4x

% Define the linear shooting function (you will need to implement this)
f = @(x, y, yp) q(x) * y + r(x);  % ODE in terms of y and y'

% Call the linear shooting function with the specified BVP for h = 1/2
[x_vals_half, y_vals_half] = linear_shooting_method(p, q, r, a, b, alpha, beta, N_half);

% Call the linear shooting function with the specified BVP for h = 1/4
[x_vals_quarter, y_vals_quarter] = linear_shooting_method(p, q, r, a, b, alpha, beta, N_quarter);

% Plot the numerical solution for h = 1/2 and h = 1/4
figure;
plot(x_vals_half, y_vals_half, 'r-o', 'DisplayName', 'Numerical Solution (h = 1/2)', 'LineWidth', 1);
hold on;
plot(x_vals_quarter, y_vals_quarter, 'b-s', 'DisplayName', 'Numerical Solution (h = 1/4)', 'LineWidth', 1);

% Define the exact solution for comparison
exact_solution = @(x) x + exp(2) * (exp(2 * x) - exp(-2 * x)) / (exp(4) - 1);

% Plot the exact solution
y_exact = arrayfun(exact_solution, x_vals_quarter);
plot(x_vals_quarter, y_exact, 'k--', 'DisplayName', 'Exact Solution', 'LineWidth', 1);

% Add labels and legend
xlabel('x');
ylabel('y');
title('Solution of BVP using Linear Shooting Method (Question2)');
legend;
grid on;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTION-03
%QUESTION-03
% Define parameters for the BVP
a = 0;  % x = 0
b = 1;  % x = 1
alpha = 1;  % y(0) = 1
beta = exp(-10);  % y(1) = e^(-10)
N_tenth = 10; % Number of steps for h = 1/10 (N = (b - a)/h = 1/10)
N_twentieth = 20; % Number of steps for h = 1/20 (N = (b - a)/h = 1/20)

% Define p, q, and r for the BVP y'' = 100y
p = @(x) 0;  % p(x) = 0
q = @(x) 100;  % q(x) = 100
r = @(x) 0;  % r(x) = 0

% Define the linear shooting function (you will need to implement this)
f = @(x, y, yp) q(x) * y + r(x);  % ODE in terms of y and y'

% Call the linear shooting function with the specified BVP for h = 1/10
[x_vals_tenth, y_vals_tenth] = linear_shooting_method(p, q, r, a, b, alpha, beta, N_tenth);

% Call the linear shooting function with the specified BVP for h = 1/20
[x_vals_twentieth, y_vals_twentieth] = linear_shooting_method(p, q, r, a, b, alpha, beta, N_twentieth);

% Plot the numerical solution for h = 1/10 and h = 1/20
figure;
plot(x_vals_tenth, y_vals_tenth, 'b-o', 'DisplayName', 'Numerical Solution (h = 1/10)', 'LineWidth', 1);
hold on;
plot(x_vals_twentieth, y_vals_twentieth, 'r-s', 'DisplayName', 'Numerical Solution (h = 1/20)', 'LineWidth', 1);

% Define the exact solution for comparison
exact_solution = @(x) exp(-10*x);

% Plot the exact solution
y_exact = arrayfun(exact_solution, x_vals_twentieth);
plot(x_vals_twentieth, y_exact, 'k--', 'DisplayName', 'Exact Solution', 'LineWidth', 1);

% Add labels and legend
xlabel('x');
ylabel('y');
title('Solution of BVP using Linear Shooting Method (Question3)');
legend;
grid on;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function for the Linear Shooting Method
function [x_vals, y_vals] = linear_shooting_method(p, q, r, a, b, alpha, beta, N)
    h = (b - a) / N; % Step size
    x_vals = linspace(a, b, N+1); % Discretized x values
    
    % Solve IVP (I): y1'' = p(x)y1' + q(x)y1 + r(x) with y1(a) = alpha, y1'(a) = 0
    [~, Y1] = runge_kutta_2nd_order(@(x, y, yp) [yp; p(x)*yp + q(x)*y + r(x)], a, [alpha; 0], h, N);
    y1_b = Y1(1, end); % y1(b)
    
    % Solve IVP (II): y2'' = p(x)y2' + q(x)y2 with y2(a) = 0, y2'(a) = 1
    [~, Y2] = runge_kutta_2nd_order(@(x, y, yp) [yp; p(x)*yp + q(x)*y], a, [0; 1], h, N);
    y2_b = Y2(1, end); % y2(b)
    
    % Combine solutions to form the final solution
    C = (beta - y1_b) / y2_b;
    y_vals = Y1(1, :) + C * Y2(1, :);
end

% Runge-Kutta solver for second-order ODE
function [x_vals, Y] = runge_kutta_2nd_order(f, x0, y0, h, N)
    x_vals = x0:h:(x0 + N*h);
    Y = zeros(2, N+1);
    Y(:, 1) = y0; % Initial values [y(0); y'(0)]
    
    for i = 1:N
        xi = x_vals(i);
        yi = Y(:, i);

        k1 = h * f(xi, yi(1), yi(2));
        k2 = h * f(xi + h/2, yi(1) + k1(1)/2, yi(2) + k1(2)/2);
        k3 = h * f(xi + h/2, yi(1) + k2(1)/2, yi(2) + k2(2)/2);
        k4 = h * f(xi + h, yi(1) + k3(1), yi(2) + k3(2));

        Y(:, i+1) = yi + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
end