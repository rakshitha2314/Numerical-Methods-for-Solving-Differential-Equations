% crank_nicolson.m
function [] = crank_nicolson(h, k, alpha, x_init, x_fin, t_init, t_req, init_cond, exact_sol)
    % Crank-Nicolson - Crank-Nicolson Method for 1D heat equation
    % h      - spatial step size
    % k      - time step size
    % alpha  - thermal diffusivity constant
    % x_init - initial x value
    % x_fin  - final x value
    % t_init - initial time
    % t_req  - required time to compute solution
    % init_cond - initial condition (u(x,0))
    % exact_sol - exact solution at the required time for comparison

    l = x_fin - x_init; % Length of the spatial domain
    x_parts = l/h; % Number of spatial partitions
    x = linspace(x_init, x_fin, x_parts+1); % Discretized x values

    t_parts = (t_req - t_init)/k; % Number of time steps

    r = k * alpha^2 / h^2; % Constant for the Crank-Nicolson method

    % Set up the coefficient matrices A and B (both tridiagonal)
    % A corresponds to the implicit part
    A = (1 + r) * eye(x_parts-1) - 0.5 * r * diag(ones(x_parts-2,1), 1) - 0.5 * r * diag(ones(x_parts-2,1), -1);

    % B corresponds to the explicit part
    B = (1 - r) * eye(x_parts-1) + 0.5 * r * diag(ones(x_parts-2,1), 1) + 0.5 * r * diag(ones(x_parts-2,1), -1);

    % Initial condition
    u = init_cond;

    % Time-stepping loop (solving the system A * u_new = B * u_old at each step)
    for n = 1:t_parts
        % Solve the system for interior points only (1D heat equation)
        u(2:x_parts) = A \ (B * u(2:x_parts)'); % Implicit solution update for interior points
    end

    % Plotting
    figure;
    grid on;
    plot(x, exact_sol, 'bo-', 'LineWidth', 2); % Plot Exact Solution in blue
    hold on;
    plot(x, u, 'ro-', 'MarkerSize', 6); % Plot Numerical Solution with black circles
    title(sprintf('Crank-Nicolson Method at t = %.2f', t_req)); % Corrected method name
    xlabel('x');
    ylabel('u(x,t)');
    legend('Exact Solution', 'Crank-Nicolson Solution'); % Corrected legend
    hold off;
end

% ====================== Main Script ======================

%QUESTION-04

%Part-(a)
% Parameters for crank_nicolson Method
alpha = 1;
x_init = 0;
x_fin = 2;
t_init = 0;
t_req = 0.5;
h = 0.4;
k = 0.1;

% Spatial domain and partitions
x_parts = (x_fin - x_init)/h;
x = linspace(x_init, x_fin, x_parts+1);

% Initial Condition (u(x, 0) = sin(pi*x))
init_cond = zeros(1, x_parts+1);
for n = 1:x_parts+1
    init_cond(n) = sin(2*pi*x(n));
end

% Exact Solution at t_req
exact_sol = zeros(1, x_parts+1);
for n = 1:x_parts+1
    exact_sol(n) = exp(-4*pi^2 * t_req) * sin(2*pi*x(n));
end

% Call the crank_nicolson function
crank_nicolson(h, k, alpha, x_init, x_fin, t_init, t_req, init_cond, exact_sol);

%Part-(b)
% Parameters for crank_nicolson Method
alpha = 1;
x_init = 0;
x_fin = 2;
t_init = 0;
t_req = 0.5;
h = 0.4;
k = 0.05;

% Spatial domain and partitions
x_parts = (x_fin - x_init)/h;
x = linspace(x_init, x_fin, x_parts+1);

% Initial Condition (u(x, 0) = sin(pi*x))
init_cond = zeros(1, x_parts+1);
for n = 1:x_parts+1
    init_cond(n) = sin(2*pi*x(n));
end

% Exact Solution at t_req
exact_sol = zeros(1, x_parts+1);
for n = 1:x_parts+1
    exact_sol(n) = exp(-4*pi^2 * t_req) * sin(2*pi*x(n));
end

% Call the crank_nicolson function
crank_nicolson(h, k, alpha, x_init, x_fin, t_init, t_req, init_cond, exact_sol);

%QUESTION-05
% Parameters for crank_nicolson Method
alpha = 1;
x_init = 0;
x_fin = pi;
t_init = 0;
t_req = 0.5;
h = pi/10;
k = 0.05;

% Spatial domain and partitions
x_parts = (x_fin - x_init)/h;
x = linspace(x_init, x_fin, x_parts+1);

% Initial Condition (u(x, 0) = sin(pi*x))
init_cond = zeros(1, x_parts+1);
for n = 1:x_parts+1
    init_cond(n) = sin(x(n));
end

% Exact Solution at t_req
exact_sol = zeros(1, x_parts+1);
for n = 1:x_parts+1
    exact_sol(n) = exp(- t_req) * sin(x(n));
end

% Call the crank_nicolson function
crank_nicolson(h, k, alpha, x_init, x_fin, t_init, t_req, init_cond, exact_sol);