% ftcs.m
function [] = ftcs(h, k, alpha, x_init, x_fin, t_init, t_req, init_cond, exact_sol)
    % FTCS - Forward Time Central Space Method for 1D heat equation
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

    u = init_cond; % Initial condition for the solution

    % Time-stepping loop
    for n = 1:t_parts
        u_new = u; % Initialize new solution array as a copy of u
        for i = 2:x_parts % Skip boundaries (1 and x_parts+1)
            u_new(i) = u(i) + k * alpha^2 / h^2 * (u(i+1) - 2*u(i) + u(i-1));
        end
        u = u_new; % Update u with the new values
    end

    % Plotting
    figure;
    grid on;
    plot(x, exact_sol, 'bo-', 'LineWidth', 0.5); % Plot Exact Solution in blue
    hold on;
    plot(x, u, 'ro-', 'MarkerSize', 6); % Plot Numerical Solution with black circles
    title(sprintf('FTCS Method at t = %.2f', t_req)); % Corrected method name
    xlabel('x');
    ylabel('u(x,t)');
    legend('Exact Solution', 'FTCS Solution'); % Corrected legend
    hold off;
end

% ====================== Main Script ======================

%QUESTION-01

% PART-(a)
% Parameters for FTCS Method
alpha = 1;
x_init = 0;
x_fin = 1;
t_init = 0;
t_req = 0.5;
h = 0.1;
k = 0.0005;

% Spatial domain and partitions
x_parts = (x_fin - x_init)/h;
x = linspace(x_init, x_fin, x_parts+1);

% Initial Condition (u(x, 0) = sin(pi*x))
init_cond = zeros(1, x_parts+1);
for n = 1:x_parts+1
    init_cond(n) = sin(pi*x(n));
end

% Exact Solution at t_req
exact_sol = zeros(1, x_parts+1);
for n = 1:x_parts+1
    exact_sol(n) = exp(-pi^2 * t_req) * sin(pi*x(n));
end

% Call the FTCS function
%ftcs(h, k, alpha, x_init, x_fin, t_init, t_req, init_cond, exact_sol);

% PART-(b)
% Parameters for FTCS Method
alpha = 1;
x_init = 0;
x_fin = 1;
t_init = 0;
t_req = 0.5;
h = 0.1;
k = 0.01;

% Spatial domain and partitions
x_parts = (x_fin - x_init)/h;
x = linspace(x_init, x_fin, x_parts+1);

% Initial Condition (u(x, 0) = sin(pi*x))
init_cond = zeros(1, x_parts+1);
for n = 1:x_parts+1
    init_cond(n) = sin(pi*x(n));
end

% Exact Solution at t_req
exact_sol = zeros(1, x_parts+1);
for n = 1:x_parts+1
    exact_sol(n) = exp(-pi^2 * t_req) * sin(pi*x(n));
end

% Call the FTCS function
%ftcs(h, k, alpha, x_init, x_fin, t_init, t_req, init_cond, exact_sol);


%QUESTION-04

%Part-(a)
% Parameters for FTCS Method
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

% Call the FTCS function
%ftcs(h, k, alpha, x_init, x_fin, t_init, t_req, init_cond, exact_sol);

%Part-(b)
% Parameters for FTCS Method
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

% Call the FTCS function
%ftcs(h, k, alpha, x_init, x_fin, t_init, t_req, init_cond, exact_sol);

%QUESTION-05
% Parameters for FTCS Method
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

% Call the FTCS function
ftcs(h, k, alpha, x_init, x_fin, t_init, t_req, init_cond, exact_sol);