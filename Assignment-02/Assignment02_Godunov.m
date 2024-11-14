% MTL 712 

% ASSIGNMENT 02 - Part01 - Godunov Method



function [] = godunov(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln, f)
    % Grid spacing
    delta_x = (x_final - x_init) / x_parts;
    delta_t = lambda * delta_x;
    t_parts = ceil((t_req - t_init) / delta_t); % Ensuring integer time steps

    % Spatial grid
    x = linspace(x_init, x_final, x_parts);
    
    % Initial conditions
    u = init_cond;
    f_u = zeros(1, x_parts); % Ensure f_u is a column vector
    for i = 1:x_parts
        f_u(i) = f(u(i));
    end

    % Time-stepping using the Godunov method
    for n = 1:t_parts
        u_new = zeros(size(init_cond));
        for i = 2:x_parts-1
            % Calculate c_plus and c_minus using the Godunov method

            if u(i) ~= u(i+1)
                % c_plus calculation
                fun_plus = @(x) (f(x) - f_u(i)) / (u(i+1) - u(i));
                mini = min(u(i), u(i+1));
                maxi = max(u(i), u(i+1));
                u_plus_star = fminbnd(fun_plus, mini, maxi);
                c_plus = -lambda * (f(u_plus_star) - f_u(i)) / (u(i+1) - u(i));
            else
                c_plus = 0; % No change if u(i) == u(i+1)
            end

            if u(i) ~= u(i-1)
                % c_minus calculation
                fun_minus = @(x) (f(x) - f_u(i)) / (u(i) - u(i-1));
                mini = min(u(i-1), u(i));
                maxi = max(u(i-1), u(i));
                u_minus_star = fminbnd(fun_minus, mini, maxi);
                c_minus = lambda * (f_u(i) - f(u_minus_star)) / (u(i) - u(i-1));
            else
                c_minus = 0; % No change if u(i) == u(i-1)
            end

            u_new(i) = u(i) + c_plus * (u(i+1) - u(i)) - c_minus * (u(i) - u(i-1));
        end

        % Periodic boundary conditions:
        if u(1) ~= u(2)
            fun_plus = @(x) (f(x) - f_u(1)) / (u(2) - u(1));
            mini = min(u(1), u(2));
            maxi = max(u(1), u(2));
            u_plus_star = fminbnd(fun_plus, mini, maxi);
            c_plus = -lambda * (f(u_plus_star) - f_u(1)) / (u(2) - u(1));
        else
            c_plus = 0;
        end

        if u(1) ~= u(x_parts)
            fun_minus = @(x) (f(x) - f_u(1)) / (u(1) - u(x_parts));
            mini = min(u(1), u(x_parts));
            maxi = max(u(1), u(x_parts));
            u_minus_star = fminbnd(fun_minus, mini, maxi);
            c_minus = lambda * (f_u(1) - f(u_minus_star)) / (u(1) - u(x_parts));
        else
            c_minus = 0;
        end

        u_new(1) = u(1) + c_plus * (u(2) - u(1)) - c_minus * (u(1) - u(x_parts));

        if u(x_parts) ~= u(1)
            fun_plus = @(x) (f(x) - f_u(x_parts)) / (u(1) - u(x_parts));
            mini = min(u(1), u(x_parts));
            maxi = max(u(1), u(x_parts));
            u_plus_star = fminbnd(fun_plus, mini, maxi);
            c_plus = -lambda * (f(u_plus_star) - f_u(x_parts)) / (u(1) - u(x_parts));
        else
            c_plus = 0;
        end

        if u(x_parts) ~= u(x_parts-1)
            fun_minus = @(x) (f(x) - f_u(x_parts)) / (u(x_parts) - u(x_parts-1));
            mini = min(u(x_parts-1), u(x_parts));
            maxi = max(u(x_parts-1), u(x_parts));
            u_minus_star = fminbnd(fun_minus, mini, maxi);
            c_minus = lambda * (f_u(x_parts) - f(u_minus_star)) / (u(x_parts) - u(x_parts-1));
        else
            c_minus = 0;
        end

        u_new(x_parts) = u(x_parts) + c_plus * (u(1) - u(x_parts)) - c_minus * (u(x_parts) - u(x_parts-1));

        % Update the solution and f_u
        u = u_new;
        for i = 1:x_parts
            f_u(i) = f(u(i));
        end
    end

    % Generate required plots:
    figure;
    plot(x, exact_soln, 'b-'); % Plot Exact Solution in blue
    hold on;
    plot(x, u, 'ko'); % Plot Numerical Solution with black circles
    title(sprintf('Godunov Method at t = %d', t_req)); % Corrected method name
    xlabel('x-val');
    ylabel('u-req');
    legend('Exact Solution', 'Godunov Solution'); % Updated legend
    
    % Optionally disable the grid
    grid off;
    
    % Optionally hide the x-axis and its ticks
    ax = gca;
    ax.XColor = 'none'; % Hide x-axis line
    ax.XTick = [];      % Remove x-axis ticks
    
    hold off;
end

% Test Case - 01

x_init = -1;
x_final = 1; 
x_parts = 40; 
lambda = 0.8; 
t_init = 0;
t_req = 30;

x = linspace(x_init, x_final, x_parts);
init_cond = -sin(pi*x); 

exact_soln = init_cond;
f = @(x) x;
 
godunov(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln,f)


% Test Case - 02

x_init = -1;
x_final = 1; 
x_parts = 40; 
lambda = 0.8; 
t_init = 0;
t_req = 4;

x = linspace(x_init, x_final, x_parts);
init_cond = zeros(1,x_parts);
for i = 1:x_parts
    if abs(x(i)) < 1/3
        init_cond(i) = 1;
    else
        init_cond(i) = 0;
    end
end

exact_soln = init_cond;
f = @(x) x;

godunov(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln, f)



% Test Case - 03

x_init = -1;
x_final = 1; 
x_parts = 600; 
lambda = 0.8; 
t_init = 0;
t_req_1 = 4;
t_req_2 = 40;

x = linspace(x_init, x_final, x_parts);
init_cond = zeros(1,x_parts);
for i = 1:x_parts
    if abs(x(i)) < 1/3
        init_cond(i) = 1;
    else
        init_cond(i) = 0;
    end
end

exact_soln = init_cond;
f = @(x) x;

godunov(x_init, x_final, x_parts, lambda, t_init, t_req_1, init_cond, exact_soln, f)
hold on;
godunov(x_init, x_final, x_parts, lambda, t_init, t_req_2, init_cond, exact_soln, f)




% Test Case - 04

x_init = -1;
x_final = 1; 
x_parts = 40; 
lambda = 0.8; 
t_init = 0;
t_req = 0.6;


x = linspace(x_init, x_final, x_parts);
init_cond = zeros(1,x_parts);
for i = 1:x_parts
    if abs(x(i)) < 1/3
        init_cond(i) = 1;
    else
        init_cond(i) = 0;
    end
end

exact_soln = zeros(1,x_parts);
for i = 1:x_parts
    if x(i) > -1/3 && x(i) < -1/3 + t_req
        exact_soln(i) = (x(i) + 1/3)/t_req;
    elseif x(i) < 1/3 + t_req/2 && x(i) > -1/3 + t_req
        exact_soln(i) = 1;
    end
end
        

f = @(x) x^2/2;

godunov(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln, f)


% Test Case - 05

x_init = -1;
x_final = 1; 
x_parts = 40; 
lambda = 0.8; 
t_init = 0;
t_req = 0.3;


x = linspace(x_init, x_final, x_parts);
init_cond = zeros(1,x_parts);
for i = 1:x_parts
    if abs(x(i)) < 1/3
        init_cond(i) = 1;
    else
        init_cond(i) = -1;
    end
end

exact_soln = zeros(1,x_parts);
for i = 1:x_parts
    if x(i) < -1/3 - t_req
        exact_soln(i) = -1;
    elseif x(i) > -1/3 - t_req && x(i) < -1/3 + t_req
        exact_soln(i) = -1 + (x(i) + 1/3 + t_req)/t_req;
    elseif x(i) < 1/3 && x(i) > -1/3 + t_req
        exact_soln(i) = 1;
    else
        exact_soln(i) = -1;
    end
end
        

f = @(x) x^2/2;

godunov(x_init, x_final, x_parts, lambda, t_init, t_req, init_cond, exact_soln, f)
