function [x_vals, y_vals] = nonlinear_shooting(f, x_range, y_boundary, slope_guess, tol, h, max_iter)
    x_start = x_range(1);
    x_end = x_range(2);
    y_start = y_boundary(1);
    y_end = y_boundary(2);
    initial_slope = slope_guess;

    for iter = 1:max_iter
        [x_vals, y] = rk4_solver(f, x_start, x_end, h, [y_start; initial_slope]);
        if abs(y(1, end) - y_end) < tol
            y_vals = y(1, :);
            return;
        end

        % Finite difference for derivative estimation
        delta = 1e-4;
        [~, y_perturbed] = rk4_solver(f, x_start, x_end, h, [y_start; initial_slope + delta]);
        derivative = (y_perturbed(1, end) - y(1, end)) / delta;
        slope_correction = (y(1, end) - y_end) / derivative;
        initial_slope = initial_slope - slope_correction;
    end
    error('Nonlinear shooting method did not converge');
end

function [x_vals, y_vals] = rk4_solver(f, x_start, x_end, h, y_init)
    x_vals = x_start:h:x_end;
    y_vals = zeros(length(y_init), length(x_vals));
    y_vals(:, 1) = y_init;

    for i = 1:length(x_vals) - 1
        x = x_vals(i);
        y = y_vals(:, i);

        k1 = h * f(x, y);
        k2 = h * f(x + h / 2, y + k1 / 2);
        k3 = h * f(x + h / 2, y + k2 / 2);
        k4 = h * f(x + h, y + k3);

        y_vals(:, i + 1) = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    end
end


%QUESTION-04

% Q4 Parameters
f = @(x, y) [y(2); (1/8)*(32 + 2*x^3 - y(1)*y(2))]; % Differential equation as a system
x_range = [1, 3];
y_boundary = [17, 43/3];
slope_guess = 0; % Initial guess for y'
tol = 1e-5;
h = (x_range(2) - x_range(1)) / 20; % Step size

% Nonlinear shooting method
[x_vals, y_vals] = nonlinear_shooting(f, x_range, y_boundary, slope_guess, tol, h, 10);

% Exact solution for Q4
y_exact = @(x) x.^2 + 16 ./ x;

% Generate exact solution values
x_exact = linspace(x_range(1), x_range(2), 100);
y_exact_vals = y_exact(x_exact);

% Plot approximate and exact solutions
figure;
plot(x_vals, y_vals, 'bo-', 'DisplayName', 'Approximate Solution', LineWidth=1);
hold on;
plot(x_exact, y_exact_vals, 'r-', 'DisplayName', 'Exact Solution', LineWidth=1);
xlabel('x');
ylabel('y');
title('Approximate vs Exact Solution (Question4)');
legend;
grid on;

%QUESTION-05

% Q5 Parameters
f = @(x, y) [y(2); -(y(2))^2 - y(1) + log(x)];
x_range = [1, 2];
y_boundary = [0, log(2)];
slope_guess = -1; % Initial guess for y'
tol = 1e-5;
h = 0.5; % Step size

% Nonlinear shooting method
[x_vals, y_vals] = nonlinear_shooting(f, x_range, y_boundary, slope_guess, tol, h, 10);

% Exact solution for Q5
y_exact = @(x) log(x);

% Generate exact solution values
x_exact = linspace(x_range(1), x_range(2), 100);
y_exact_vals = y_exact(x_exact);

% Plot approximate and exact solutions
figure;
plot(x_vals, y_vals, 'bo-', 'DisplayName', 'Approximate Solution', LineWidth=1);
hold on;
plot(x_exact, y_exact_vals, 'r-', 'DisplayName', 'Exact Solution', LineWidth=1);
xlabel('x');
ylabel('y');
title('Approximate vs Exact Solution (Question5)');
legend;
grid on;

%QUESTION-06

% Q6 Parameters
f = @(x, y) [y(2); -exp(-2*y(1))];
x_range = [1, 2];
y_boundary = [0, log(2)];
slope_guess = -1; % Initial guess for y'
tol = 1e-4;
h = (x_range(2) - x_range(1)) / 10; % Step size

% Nonlinear shooting method
[x_vals, y_vals] = nonlinear_shooting(f, x_range, y_boundary, slope_guess, tol, h, 10);

% Exact solution for Q6
y_exact = @(x) log(x);

% Generate exact solution values
x_exact = linspace(x_range(1), x_range(2), 100);
y_exact_vals = y_exact(x_exact);

% Plot approximate and exact solutions
figure;
plot(x_vals, y_vals, 'bo-', 'DisplayName', 'Approximate Solution', LineWidth=1);
hold on;
plot(x_exact, y_exact_vals, 'r-', 'DisplayName', 'Exact Solution', LineWidth=1);
xlabel('x');
ylabel('y');
title('Approximate vs Exact Solution (Question6)');
legend;
grid on;

% QUESTION-07
% Q7 Parameters
f = @(x, y) [y(2); y(2)*cos(x) - y(1)*log(y(1))];
x_range = [0, pi/2];
y_boundary = [1, exp(1)];
slope_guess = 1; % Initial guess for y'
tol = 1e-5;
h = (x_range(2) - x_range(1)) / 10; % Step size

% Nonlinear shooting method
[x_vals, y_vals] = nonlinear_shooting(f, x_range, y_boundary, slope_guess, tol, h, 10);

% Exact solution for Q7
y_exact = @(x) exp(sin(x));

% Generate exact solution values
x_exact = linspace(x_range(1), x_range(2), 100);
y_exact_vals = y_exact(x_exact);

% Plot approximate and exact solutions
figure;
plot(x_vals, y_vals, 'bo-', 'DisplayName', 'Approximate Solution', LineWidth=1);
hold on;
plot(x_exact, y_exact_vals, 'r-', 'DisplayName', 'Exact Solution', LineWidth=1);
xlabel('x');
ylabel('y');
title('Approximate vs Exact Solution (Question7)');
legend;
grid on;